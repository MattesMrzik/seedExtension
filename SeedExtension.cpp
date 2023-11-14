#include "SeedExtension.hpp"

#include "MetagraphInterface.h"
#include "IdentifierMapping.h"
#include "Link.h"
#include "PathBundleTip.hpp"

#include <iostream>
#include <string>
#include <vector>

void SeedExtension::initFirstTip(std::vector<MetagraphInterface::NodeID> nodeIDs,
                                 LinkPtr link,
                                 std::shared_ptr<IdentifierMapping const> idMap_) {
    idMap = idMap_; //not in ctor, bc idMap not available in test/testSeedExtension.cpp
    tipsHistory.clear();
    upStreamTipsHistory.clear();

    // for extension analysis
    tooManyDeletedAnnos = 0;
    nMerges = 0;
    nSplits = 0;
    tooManyAnnosInInit = 0;
    maxUpstreamSteps = 0;
    maxSteps = 0;

    auto firstAllTips = std::make_shared<AllTips>(AllTips{});

    //bc nodeIDs contains ids twice, see implementation in Linkset.h
    std::unordered_set<MetagraphInterface::NodeID> nonDupNodeIDs;
    nonDupNodeIDs.insert(nodeIDs.begin(), nodeIDs.end());

    for (auto nodeID : nonDupNodeIDs) {
        auto firstTip = std::make_shared<PathBundleTip>(PathBundleTip{});
        firstTip->nodeID = nodeID;

        //only consider annos in graph which were provided by link = seed
        for (auto & occurrence : link->occurrence()) {
            for (auto & metaAnno : graph->getAnnotation(nodeID)) {
                if (idMap->queryGenomeName(occurrence.genome()) == metaAnno.genome &&
                    idMap->querySequenceName(occurrence.sequence()) == metaAnno.sequence &&
                    occurrence.reverse() == metaAnno.reverse_strand &&
                    occurrence.position() == metaAnno.bin_idx) {

                    firstTip->annotations.insert({metaAnno,{0,0,0,0}});
                }
            }
        }
        firstAllTips->tips.insert({nodeID, firstTip});
    }

    firstAllTips->numberOfExtensionsMade = 0;
    firstAllTips->totalScore = 0;

    tooManyAnnosInInit = firstAllTips->nAnnotations() - link->occurrence().size();

    tipsHistory.push_back(firstAllTips);
    upStreamTipsHistory.push_back(std::make_shared<AllTips>(*firstAllTips));
    // score of initial kmners assigned to downStream annos
    // dont initScore for upStream, bc then score of initial kmer would count twice to totalScore
    tipsHistory.back()->initScore(graph);
}
void SeedExtension::extend(size_t sufficientMaxScore) {
    // extend downStream
    extendOneSide(sufficientMaxScore, tipsHistory, false);
    // extend upStream
    extendOneSide(sufficientMaxScore, upStreamTipsHistory, true);

}

void SeedExtension::extendOneSide(size_t sufficientMaxScore,
                                  std::vector<std::shared_ptr<AllTips>> & tipsHis,
                                  bool upStream) {
    while (tipsHis.back()->nGenomes() >= 2 &&
           totalScore()  < (int)sufficientMaxScore &&
           tipsHis.back()->containsReferenzGenome(idMap) &&
           tipsHis.size() < sufficientMaxScore * 3) { //catches potential inf loop

        std::vector<int> splitsAndMerge{0,0}; // collects some info about extension
        auto newStep = tipsHis.back()->extendAllTipsWithAnalysis(graph, upStream, binsize, splitsAndMerge);
        nSplits += splitsAndMerge[0];
        nMerges += splitsAndMerge[1];

        tipsHis.push_back(newStep); // add new AllTips to back of Alignment
        if (upStream) {
            maxUpstreamSteps = newStep->numberOfExtensionsMade;
        }
        else {
            maxSteps = newStep->numberOfExtensionsMade;
        }

        // delete some annos if necessary
        xDrop(config->xdrop(), tipsHis, upStream);
    }
}
//! find the annos in current allTips which have a score less than:
//! their their maxscore - xdrop
//! and are furthest in the past
size_t SeedExtension::getAnnosToBeDropped(uint64_t xdrop,
                                     std::shared_ptr<AllTips> allTips,
                                     std::vector<MetagraphInterface::NodeAnnotation> & annosToBeDropped) const {
    auto xDropFurthestBack = allTips->numberOfExtensionsMade;
    for (auto & [nodeID, tip] : allTips->tips) {
        for (auto & [metaAnno, annoScore] : tip->annotations) {
            if (annoScore.currentScore < annoScore.maxScore - xdrop ) {
                // delete only the xdrop with the maxscore furthest in the past
                if (annosToBeDropped.size() == 0 ||
                    annoScore.ageOfMaxScore == xDropFurthestBack) {

                    annosToBeDropped.push_back(metaAnno);
                    xDropFurthestBack = annoScore.ageOfMaxScore; // for the case if annosToBeDropped.size() == 0

                }
                // only drop annotations which have the oldest (=lowest) ageOfMaxScore
                else if (annoScore.ageOfMaxScore < xDropFurthestBack) {
                    annosToBeDropped.clear();
                    annosToBeDropped.push_back(metaAnno);
                    xDropFurthestBack = annoScore.ageOfMaxScore;
                }
            }
        }
    }
    return(xDropFurthestBack);
}
// cant just do multi pop here, since when xdrop triggers, a copy of allTips is made
// and pushed to vector but with identical numberOfExtensionsMade
// therefore vector index is not always equal to vector[index].numberOfExtensionsMade
// if annos are deleted in place and no copy of alltips is made
// then multi pop could be used instead of this funtion
void SeedExtension::undoSteps(std::vector<std::shared_ptr<AllTips>> & tipsHis,
                              size_t goBackTo) {
    uint64_t latestNumberOfExtensionsMade = tipsHis.back()->numberOfExtensionsMade;
    while (latestNumberOfExtensionsMade != goBackTo) {
        tipsHis.pop_back();
        latestNumberOfExtensionsMade = tipsHis.back()->numberOfExtensionsMade;
    }
}
std::vector<unsigned>
SeedExtension::removeAnnos(tipsMapType & tips,
                           std::vector<MetagraphInterface::NodeAnnotation> & annosToBeDropped,
                           bool upStream,
                           size_t nBackSteps) {
    // to update score after annos have been removed
    std::vector<unsigned> acgt{0,0,0,0}; // will be returned

    for (auto const & annoToBeDropped : annosToBeDropped) {
        for (auto tipItr = tips.begin(); tipItr != tips.end(); ) {
            auto & tip = tipItr->second;
            char base = upStream ?
                        graph->getKmer(tip->nodeID).front() :
                        graph->getKmer(tip->nodeID).back();
            auto foundAnnoPtr = tip->annotations.find(annoToBeDropped);
            if (foundAnnoPtr != tip->annotations.end()) {
                //add base to acgt to return, to update score after xdrop
                acgt[AllTips::baseToId(base)] += 1;
                // found annoToBeDropped in tip
                tip->annotations.erase(foundAnnoPtr);
            }
            else {
                // didnt find annoToBeDropped in tip
                // -> look for annoToBeDropped with modified bin_idx
                const int signedBinsize = upStream ? -binsize : binsize;
                // how many transition could have been made in nBackSteps extension steps
                unsigned maxNumberOfBinSizeTransitions = (nBackSteps - 1) / binsize + 1;
                unsigned minNumberOfBinSizeTransitions = nBackSteps / binsize;
                for (unsigned binsizesback = minNumberOfBinSizeTransitions;
                        binsizesback <= maxNumberOfBinSizeTransitions;
                        binsizesback++) { // this loop is at min 1 and at max 2 cycles long
                    // check whether modified bin_idx is legal
                    int tilediff = binsizesback * signedBinsize;
                    if (tilediff > 0 && annoToBeDropped.bin_idx < (unsigned)tilediff) {
                        // not legal -> skip it
                        continue;
                    }
                    // modified annoToBeDropped used for searching anno with different bin_idx
                    MetagraphInterface::NodeAnnotation upStreamAnno {annoToBeDropped.genome,
                                                                     annoToBeDropped.sequence,
                                                                     annoToBeDropped.reverse_strand,
                                                                     annoToBeDropped.bin_idx - tilediff};
                    auto foundUpStreamAnnoPtr = tip->annotations.find(upStreamAnno);
                    if (foundUpStreamAnnoPtr != tip->annotations.end()) {
                        acgt[AllTips::baseToId(base)] += 1;
                        // found modified annoToBeDropped -> delete it
                        tip->annotations.erase(foundUpStreamAnnoPtr);
                    }
                }
            }
            // if tip is empty -> delete it
            if (tip->annotations.size() == 0) {
                tipItr = tips.erase(tipItr);
            }
            else {
                tipItr++;
            }
        }
    }
    return acgt;
}
void SeedExtension::xDrop(uint64_t xdrop,
                          std::vector<std::shared_ptr<AllTips>> & tipsHis,
                          bool upStream) {

    std::vector<MetagraphInterface::NodeAnnotation> annosToBeDropped;
    auto goBackTo = getAnnosToBeDropped(xdrop, tipsHis.back(), annosToBeDropped);
    // if goBackTo is close to nodes where extension started
    // go back to start instead
    goBackTo = goBackTo < 5 ? 0 : goBackTo;

    size_t nBackSteps = tipsHis.back()->numberOfExtensionsMade - goBackTo;

    if (annosToBeDropped.size() != 0) {

        undoSteps(tipsHis, goBackTo);

        auto nAnnotationsBeforeXDrop = tipsHis.back()->nAnnotations();

        //not trimmed yet, but when removeAnnos is called
        auto trimmedAllTips = std::make_shared<AllTips>(*tipsHis.back());
        unsigned sizeBeforeXdrop = trimmedAllTips->nAnnotations();

        auto nACGTBeforeRemove = trimmedAllTips->nACGT(upStream, graph);

        auto removedChars = removeAnnos(trimmedAllTips->tips, annosToBeDropped, upStream, nBackSteps);

        auto nACGTAfterRemove = trimmedAllTips->nACGT(upStream, graph);

        // for extension analysis
        if (sizeBeforeXdrop >= trimmedAllTips->nAnnotations() + annosToBeDropped.size()) {
            tooManyDeletedAnnos += sizeBeforeXdrop - trimmedAllTips->nAnnotations() - annosToBeDropped.size();
        }

        //update score after annos were removed
        for (auto & [id, tip] : trimmedAllTips->tips) {
            char currentBase = upStream ?
                               graph->getKmer(id).front() :
                               graph->getKmer(id).back();
            for (auto & [anno, scoreStructure] : tip->annotations) {
                auto score = trimmedAllTips->charVsProfileScore(currentBase, nACGTBeforeRemove);
                scoreStructure.currentScore -= score;
                score = trimmedAllTips->charVsProfileScore(currentBase, nACGTAfterRemove);
                scoreStructure.currentScore += score;
            }
        }

        tipsHis.push_back(trimmedAllTips);

        if (nAnnotationsBeforeXDrop != nACGTBeforeRemove[0] + nACGTBeforeRemove[1] + nACGTBeforeRemove[2] + nACGTBeforeRemove[3]) {
            std::cout << "something went wrong SeedExtension nAnnotationsBeforeXDrop" << '\n';
            exit(1);
        }
        auto nAnnotationsAfterXDrop = tipsHis.back()->nAnnotations();
        if (nAnnotationsAfterXDrop != nACGTAfterRemove[0] + nACGTAfterRemove[1] + nACGTAfterRemove[2] + nACGTAfterRemove[3]) {
            std::cout << "something went wrong SeedExtension nAnnotationsAfterXDrop" << '\n';
        }
    }
}
