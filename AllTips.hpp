#ifndef _ALLTIPS_HPP_
#define _ALLTIPS_HPP_

#include "PathBundleTip.hpp"
#include "MetagraphInterface.h"
#include "IdentifierMapping.h"

#include <vector>
#include <memory>

/* ! Represents all nodes with respective bundles of annotations
* (Erweiterungsstand)
*/
class AllTips {
public:
    using tipsMapType = std::unordered_map<uint64_t, std::shared_ptr<PathBundleTip>>;

    AllTips(AllTips const & other):numberOfExtensionsMade{other.numberOfExtensionsMade},
                                   totalScore{other.totalScore},
                                   tips{} {
        for(auto idTipPair : other.tips) {
            tips.insert({idTipPair.first, std::make_shared<PathBundleTip>(*idTipPair.second)});
        }
    }
    AllTips():numberOfExtensionsMade{},totalScore{},tips{} {}

    AllTips(uint64_t numberOfExtensionsMade_, double totalScore_):
            numberOfExtensionsMade{numberOfExtensionsMade_},
            totalScore{totalScore_},tips{} {}

    AllTips(uint64_t numberOfExtensionsMade_, double totalScore_, tipsMapType & tips_):
            numberOfExtensionsMade{numberOfExtensionsMade_},
            totalScore{totalScore_},tips{tips_} {}

    //! extends all the PathBundleTip s
    std::shared_ptr<AllTips> extendAllTips(std::shared_ptr<MetagraphInterface const> graph,
                                           bool upStream,
                                           size_t binsize);
    //! same as extendAllTips except it collects some statistics about extension
    std::shared_ptr<AllTips> extendAllTipsWithAnalysis(std::shared_ptr<MetagraphInterface const> graph,
                                                       bool upStream,
                                                       size_t binsize,
                                                       std::vector<int> & splitsAndMerge);
    //! exntends all PathBundleTip s without updating the score
    /*! is called by extendAllTips()
     * the initial scores of extended tips are idetical to current scores
     * they are then updated by updateScores()
     */
    void extendWithoutUpdatingScore(std::shared_ptr<AllTips> newAllTips,
                                    bool upStream,
                                    size_t binsize,
                                    std::shared_ptr<MetagraphInterface const>  graph) const;
    void extendWithoutUpdatingScoreWithAnalysis(std::shared_ptr<AllTips> newAllTips,
                                    bool upStream,
                                    size_t binsize,
                                    std::shared_ptr<MetagraphInterface const>  graph,
                                    std::vector<int> & splitsAndMerge) const;
    //! updates the scores of all annos and totalScore after extension without updating score
    void updateScores(bool upStream,
                      std::shared_ptr<MetagraphInterface const>  graph,
                      double previousTotalScore);

    //TODO make this static when scoring matrix is static
    static double charVsProfileScore(char b, std::vector<unsigned> & acgt);

    void initScore(std::shared_ptr<MetagraphInterface const>  graph);

    //! returns number of each base
    /*! at first or last position of all kmers corresponding to all PathBundleTip s nodeID
     * depending on extension direction
     */
    std::vector<unsigned> nACGT(bool upStream,
                            std::shared_ptr<MetagraphInterface const>  graph) const;

    //! returns number of each base
    /*! at position of all kmers corresponding to all PathBundleTip s nodeID */
    std::vector<unsigned>
    nACGTatKmersPos(unsigned pos,
                    std::shared_ptr<MetagraphInterface const>  graph) const;

    void printAllTips(std::shared_ptr<MetagraphInterface const> graph) const;

    void printAllTips() const;

    static int getScore(char base1, char base2) {
        return(scoringMatrix[baseToId(base1)][baseToId(base2)]);
    }

    bool containsReferenzGenome(std::shared_ptr<IdentifierMapping const> idMap) const;

    bool containsGenome(unsigned genomeID, std::shared_ptr<IdentifierMapping const> idMap) const;

    static unsigned baseToId(char const base);

    //! returns the sum of all annotations of all the PathBundleTip s of this AllTips
    size_t nAnnotations() const;
    //! returns the number of unique genomes
    unsigned nGenomes() const;
    //TODO maybe also for nGenomes and nSeq

    //! number of extension steps made in this direction since initFirstTip
    uint64_t numberOfExtensionsMade;
    //! the sum of scores of all annotations of all PathBundleTip s in this AllTips
    double totalScore;
    //! all the PathBundleTips
    tipsMapType tips;

    static std::vector<std::vector<int>> scoringMatrix;

};
#endif //_ALLTIPS_HPP_
