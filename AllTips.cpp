#include "AllTips.hpp"

#include "PathBundleTip.hpp"
#include "MetagraphInterface.h"

#include <iostream>
#include <vector>
#include <memory>
#include <math.h>

std::vector<std::vector<int>> AllTips::scoringMatrix{{5,-4,-4,-4},
                                                     {-4,5,-4,-4},
                                                     {-4,-4,5,-4},
                                                     {-4,-4,-4,5}};

shared_ptr<AllTips>
AllTips::extendAllTips(std::shared_ptr<MetagraphInterface const>  graph,
                               bool upStream,
                               size_t binsize) {
    std::vector<int> v {0,0};
    return(extendAllTipsWithAnalysis(graph, upStream, binsize, v));
}

shared_ptr<AllTips>
AllTips::extendAllTipsWithAnalysis(std::shared_ptr<MetagraphInterface const>  graph,
                               bool upStream,
                               size_t binsize,
                               std::vector<int> & splitsAndMerge) {
    auto newAllTips = std::make_shared<AllTips>(AllTips{});

    // extend every tip and match every annotation and then merge all tips with same id
    extendWithoutUpdatingScoreWithAnalysis(newAllTips, upStream, binsize, graph, splitsAndMerge);
    newAllTips->updateScores(upStream, graph, totalScore);

    return newAllTips;
}
void AllTips::extendWithoutUpdatingScore(std::shared_ptr<AllTips> newAllTips,
                                         bool upStream,
                                         size_t binsize,
                                         std::shared_ptr<MetagraphInterface const>  graph) const {
    std::vector<int> v {0,0};
    extendWithoutUpdatingScoreWithAnalysis(newAllTips, upStream, binsize, graph, v);
}

void AllTips::extendWithoutUpdatingScoreWithAnalysis(std::shared_ptr<AllTips> newAllTips,
                                         bool upStream,
                                         size_t binsize,
                                         std::shared_ptr<MetagraphInterface const>  graph,
                                         std::vector<int> & splitsAndMerge) const {

    for (auto & [id, tip] : tips) {
        // vector<std::shared_ptr<PathBundleTip>>
        auto const newTips = tip->extendTip(graph, upStream, binsize, numberOfExtensionsMade);
        // at max 4
        for (auto & newTip : newTips) {
            // if that node is already in allNewTips -> merging the annotations
            auto foundSameNodeID = newAllTips->tips.find(newTip->nodeID);
            if (foundSameNodeID != newAllTips->tips.end()) {
                unsigned tip1size = newTip->annotations.size();
                unsigned tip2size = foundSameNodeID->second->annotations.size();
                // merge annotations:
                for (auto & annotation : newTip->annotations) {
                    //TODO what if intersection of annotations of both tips is not empty?
                    foundSameNodeID->second->annotations.insert(annotation);
                }
                splitsAndMerge[1] = tip1size + tip2size - foundSameNodeID->second->annotations.size(); // merge
            }
            // that node isnt in allNewTips -> add it to newTips
            else {
                newAllTips->tips.insert({newTip->nodeID, newTip});
            }
        }
    }
    unsigned nAnnosBefore = nAnnotations();
    unsigned nAnnosAfter = 0;
    for (auto & [id, tip] : newAllTips->tips) {
        nAnnosAfter += tip->annotations.size();
    }
    if (nAnnosAfter >= nAnnosBefore + splitsAndMerge[1] ) {
        splitsAndMerge[0] = nAnnosAfter - nAnnosBefore + splitsAndMerge[1]; // splits
        // if (splitsAndMerge[0] + splitsAndMerge[1] != 0) {
        //     std::cout << "splits (" << splitsAndMerge[0] << "), merge[1] (" << splitsAndMerge[1] << ") " << '\n';
        // }
    }
    newAllTips->numberOfExtensionsMade = numberOfExtensionsMade + 1;
}
//! for the first AllTips evaluate every base of kmers corresponding to tips
void AllTips::initScore(std::shared_ptr<MetagraphInterface const>  graph) {
    for (unsigned i = 0; i < graph->getK(); i++) {
        auto acgt = nACGTatKmersPos(i, graph);
        for (auto & [id, tip] : this->tips) {
            auto currentBase = graph->getKmer(id).at(i);
            for (auto & [anno, scoreStruct] : tip->annotations) {
                auto score = charVsProfileScore(currentBase, acgt);
                scoreStruct.currentScore += score;
                scoreStruct.maxScore += score;
                totalScore += score;
            }
        }
    }
}

//! returns the score for base b with regards to profile
double AllTips::charVsProfileScore(char base, std::vector<unsigned> & acgt) {
    unsigned nAnnotations = acgt[0] + acgt[1] + acgt[2] + acgt[3];

    // -1 bc we dont want to score base to itself
    if (acgt[baseToId(base)] == 0) {
        std::cout << "illegal base (" << base << ") in updateAnnoScore for array acgt" << '\n';
    }
    acgt[baseToId(base)] -= 1;
    double score = scoringMatrix[baseToId(base)][baseToId('A')] * (double)(acgt[0])
                 + scoringMatrix[baseToId(base)][baseToId('C')] * (double)(acgt[1])
                 + scoringMatrix[baseToId(base)][baseToId('G')] * (double)(acgt[2])
                 + scoringMatrix[baseToId(base)][baseToId('T')] * (double)(acgt[3]);
    //undo modification
    acgt[baseToId(base)] += 1;
    return score / (double)nAnnotations;
}
/*! for one extension step update the score of all annos in allNewTips
*/
void AllTips::updateScores(bool upStream,
                           std::shared_ptr<MetagraphInterface const>  graph,
                           double previousTotalScore){
    std::vector<unsigned> acgt = nACGT(upStream, graph);
    // the delta of totalScore when extending
    double deltaScore = 0;
    for (auto & [nodeID, tip] : tips){
        char currentBase = upStream ?
                           graph->getKmer(nodeID).front() :
                           graph->getKmer(nodeID).back();
        for (auto & [metaAnno, annoScore] : tip->annotations) {
            double oldScore = annoScore.currentScore;
            double newScore = oldScore + charVsProfileScore(currentBase, acgt);
            // to update totalScore
            deltaScore += newScore - oldScore;
            // set new Score
            annoScore.currentScore = newScore;
            // update maxScore if necessary
            if (annoScore.maxScore <= annoScore.currentScore) {
                annoScore.maxScore = annoScore.currentScore;
                // set age of maxScore
                annoScore.ageOfMaxScore = numberOfExtensionsMade;
            }
        }
    }
    totalScore = previousTotalScore + deltaScore;
}
/*! returns the number of each base at position most left (= 0) if upStream and most right (= k - 1) if downStream
* of all kmers corresponding to all tips in this allNewTips
*/
std::vector<unsigned>
AllTips::nACGT(bool upStream, std::shared_ptr<MetagraphInterface const>  graph) const {
    return(nACGTatKmersPos(upStream ? 0 : (graph->getK() - 1) , graph));
}
//! returns the number of each base at position pos of all kmers corresponding to all tips in this allNewTips
std::vector<unsigned>
AllTips::nACGTatKmersPos(unsigned pos,
                         std::shared_ptr<MetagraphInterface const>  graph) const {
    if (pos >= graph->getK()) {
        std::cout << "kmer out of bounds in AllTips::nACGTatKmersPos()" << '\n';
        exit(1);
    }
    // std::cout << "nACGTatKmersPos" << '\n';
    std::vector<unsigned> acgt{0,0,0,0};
    for (auto & [nodeID, tip] : tips) {
        char base = graph->getKmer(nodeID).at(pos);
        acgt[AllTips::baseToId(base)] += tip->annotations.size();
        // if (nodeID == 10433) std::cout << "base = " << base << ", tip->annotations.size() = " << tip->annotations.size() << '\n';
    }
    return(acgt);
}
// returns the number of annotations
// |(Annos(T,e))|
size_t AllTips::nAnnotations() const {
    size_t nAnnos = 0;
    for (auto & [id, tip] : tips) {
        nAnnos += tip->annotations.size();
    }
    return(nAnnos);
}
// returns the number of unique genomes
unsigned AllTips::nGenomes() const {
    std::unordered_set<std::string> genomes;
    for (auto & [id, tip] : tips) {
        for (auto & [anno, score] : tip->annotations) {
            genomes.insert(anno.genome);
        }
    }
    return genomes.size();
}
// returns true, if config()->genome1() is present
bool AllTips::containsReferenzGenome(std::shared_ptr<IdentifierMapping const> idMap) const {
    return(containsGenome(0, idMap));
}

// returns true, genome with id <genomeID> is present
bool AllTips::containsGenome(unsigned genomeID, std::shared_ptr<IdentifierMapping const> idMap) const {
    for (auto & [id, tip] : tips) {
        for (auto & [anno, score] : tip->annotations) {
            if (anno.genome == idMap->queryGenomeName(genomeID)) {
                return true;
            }
        }
    }
    return false;
}

void AllTips::printAllTips(std::shared_ptr<MetagraphInterface const> graph) const {
    std::cout<<"==========> printing AllTips.cpp <=========="<<std::endl;
    std::cout << "numberOfExtensionsMade: " << numberOfExtensionsMade<<std::endl;
    auto numNodes = graph->numNodes();
    auto maxDecimalPlaces = std::to_string(numNodes).size();
    for (auto & [nodeID, tip] : tips) {
        auto currentDecimalPlaces = std::to_string(nodeID).size();
        int filler = maxDecimalPlaces < currentDecimalPlaces ? 0 : maxDecimalPlaces - currentDecimalPlaces;
        for (auto & [metaAnno, annoScore] : tip->annotations) {
            std::cout << nodeID;
            for (int i = 0; i < filler + 1; i++) {
                std::cout << " ";
            }
            std::cout << "label: " << graph->getKmer(nodeID)
                      << ",\tannotation: " << metaAnno.genome
                      << ", cords: " << metaAnno.bin_idx << ", ";
            int filler = std::to_string(metaAnno.bin_idx).size() > 7 ? 0 : 7 - std::to_string(metaAnno.bin_idx).size();
            for (int i = 0; i < filler; i++) {
                std::cout << " ";
            }
            std::cout << "score: " << roundf(annoScore.currentScore*100)/100
                      << ",\tmaxscore: " << roundf(annoScore.maxScore*100)/100
                      << ",\tageOfMaxScore: " << annoScore.ageOfMaxScore
                      << ",\tsequence: " <<  metaAnno.sequence
                      << ",\tnAnnotations(): " << nAnnotations()
                      << std::endl;
        }
    }
    std::cout<<"==========> printing AllTips.cpp done <====="<<std::endl;
}

void AllTips::printAllTips() const {
    std::cout<<"==========> printing AllTips.cpp <=========="<<std::endl;
    std::cout << "numberOfExtensionsMade: " << numberOfExtensionsMade<<std::endl;
    for (auto & [nodeID, tip] : tips) {
        for (auto & [metaAnno, annoScore] : tip->annotations) {
            std::cout << "NodeID:\t" << nodeID
                      << ",\tannotation: " << metaAnno.genome
                      << ", cords: " << metaAnno.bin_idx
                      << ",\tscore: " << roundf(annoScore.currentScore*100)/100
                      << ",\tmaxscore: " <<  roundf(annoScore.maxScore*100)/100
                      << std::endl;
        }
    }
    std::cout<<"==========> printing AllTips.cpp done <====="<<std::endl;
}
// converts base to id
unsigned AllTips::baseToId(char const base) {
    switch (base) {
        case 'A':return 0;
        break;
        case 'C':return 1;
        break;
        case 'G':return 2;
        break;
        case 'T':return 3;
        break;
        case 'a':return 0;
        break;
        case 'c':return 1;
        break;
        case 'g':return 2;
        break;
        case 't':return 3;
        break;
    }
    return(0);
}
