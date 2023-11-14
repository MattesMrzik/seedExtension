#include "PathBundleTip.hpp"

#include "MetagraphInterface.h"

#include <vector>
#include <memory>

//! extends every annotation in current bundle
std::vector<std::shared_ptr<PathBundleTip>>
PathBundleTip::extendTip(std::shared_ptr <MetagraphInterface const> graph,
                         bool upStream,
                         size_t binsize_,
                         uint64_t numberOfExtensionsMade) {
    const int binsize = upStream ? - binsize_ : binsize_;

    // bundles of each outgoing node, which has at least one continuing  annotation
    std::vector<std::shared_ptr<PathBundleTip>> outgoingTips;

    std::vector<long unsigned int> outgoing_ids;
    if (!upStream) outgoing_ids = graph->getOutgoing(nodeID);
    else outgoing_ids = graph->getIncoming(nodeID);

    // loop through new nodes: max 4 different
    for (auto out_id : outgoing_ids) {
        auto outgoingNodeAnnotations = graph->getAnnotation(out_id);
        auto outgoingTip = std::make_shared<PathBundleTip>(PathBundleTip{});
        // loop through annotations of new node
        for (MetagraphInterface::NodeAnnotation & outgoingNodeAnnotation : outgoingNodeAnnotations) {
            // look for the outgoingNodeAnnotation in the current PathBundleTip
            auto matchingAnno = annotations.find(outgoingNodeAnnotation);
            unsigned latestTransition;

            if (matchingAnno != annotations.end()) {
                latestTransition = matchingAnno->second.latestTransition;
                outgoingTip->annotations.insert({outgoingNodeAnnotation,
                                              {matchingAnno->second.currentScore,
                                               matchingAnno->second.maxScore,
                                               matchingAnno->second.ageOfMaxScore,
                                               latestTransition}});
            }
            // check for anno in neighbouring bin_idx
            // check whether neighbouring bin_idx isnt out of bounds
            bool modifiedBinIdxIsGreater0 = upStream ?
                                            true :
                                            binsize_ <= outgoingNodeAnnotation.bin_idx;

            if (modifiedBinIdxIsGreater0) {
                // look for outgoingNodeAnnotation in current tip
                // with modified bin_idx
                outgoingNodeAnnotation.bin_idx -= binsize;
                matchingAnno = annotations.find(outgoingNodeAnnotation);
                outgoingNodeAnnotation.bin_idx += binsize; //undo modification

                if (matchingAnno != annotations.end()) {
                    // found annotation in neighbouring bin_idx
                    // therefore a transition from anno.position to anno-position + binsize was made
                    latestTransition = numberOfExtensionsMade + 1;
                    // if the previous transition isnt old enough
                    // then this transition isnt plausible with linear sequence
                    bool latestTransitionIsOldEnough = matchingAnno->second.latestTransition == 0
                                                     || latestTransition
                                                     - matchingAnno->second.latestTransition
                                                     >= binsize_;
                    if (latestTransitionIsOldEnough) {
                        outgoingTip->annotations.insert({outgoingNodeAnnotation,
                                                      {matchingAnno->second.currentScore,
                                                       matchingAnno->second.maxScore,
                                                       matchingAnno->second.ageOfMaxScore,
                                                       latestTransition}});
                    }
                }
            }
        }
        // if no continuing annotations found -> no need to create
        // new PathBundleTip for that node
        if (outgoingTip->annotations.size() > 0){
            outgoingTip->nodeID = out_id;
            outgoingTips.push_back(outgoingTip);
        }
    }
    return outgoingTips;
}

void  PathBundleTip::print(std::shared_ptr<MetagraphInterface const> graph) const {
    for (auto & [metaAnno, annoScore] : annotations){
        std::cout<< graph->getKmer(nodeID)
                 << metaAnno.genome << " "
                 << metaAnno.bin_idx << " "
                 << annoScore.currentScore << " "
                 << annoScore.maxScore << std::endl;
    }
}

void  PathBundleTip::print() const{
    for (auto & [metaAnno, annoScore] : annotations){
        std::cout << metaAnno.genome <<" "
                  << metaAnno.bin_idx <<" "
                  << annoScore.currentScore <<" "
                  << annoScore.maxScore << std::endl;
    }
}
void PathBundleTip::printMetaAnno(MetagraphInterface::NodeAnnotation const & anno) {
    std::cout << anno.genome << ", "
              << anno.sequence << ", "
              << anno.reverse_strand << ", "
              << anno.bin_idx <<  '\n';
}
