#ifndef _SeedExtension_HPP_
#define _SeedExtension_HPP_

#include "Configuration.h"
#include "IdentifierMapping.h"
#include "MetagraphInterface.h"
#include "Link.h"
#include "AllTips.hpp"


#include <iostream>
#include <memory>
#include <vector>

/* Holds a vector of AllTips for upStream and downStream
* Represents an alignment in genom grapg
* ((T_0,e_0), (T_1,e_1), ... , (T_n,e_n))
*/
class SeedExtension{

public:
    using tipsMapType = std::unordered_map<uint64_t, std::shared_ptr<PathBundleTip>>;

    SeedExtension(size_t binsize_):tipsHistory{},
                    upStreamTipsHistory{},
                    graph{},
                    config{},
                    idMap{},
                    binsize{binsize_} {};

    SeedExtension(std::shared_ptr<MetagraphInterface const> graph_,
                  std::shared_ptr<Configuration const> config_,
                  size_t binsize_):
               tipsHistory{},
               upStreamTipsHistory{},
               graph{graph_},
               config{config_},
               idMap{},
               binsize{binsize_}{}

    //! calls the extention to both sides (upstream and downstream)
    void extend(size_t sufficientMaxScore);
    void extendOneSide(size_t sufficientMaxScore,
                       std::vector<std::shared_ptr<AllTips>> & tipsHis,
                       bool upStream);
    //! for a given seed = link, init the first Alltips (T_0,e_0)
    void initFirstTip(std::vector<MetagraphInterface::NodeID> nodeIDs,
                      LinkPtr link,
                      std::shared_ptr<IdentifierMapping const> idMap);
    //! removes the annos provided by annosToBeDropped from allTips
    std::vector<unsigned>
    removeAnnos(tipsMapType & allTips,
                std::vector<MetagraphInterface::NodeAnnotation> & annosToBeDropped,
                bool upStream,
                size_t nBackSteps);

    //! undoes the last extension steps until the goBackTo extension step is reached
    // goBackTo = a.m^\star
    void undoSteps(std::vector<std::shared_ptr<AllTips>> & tipsHis,
                                  size_t goBackTo);

    //! determine the annos whos score is less than their max score - xDrop
    size_t
    getAnnosToBeDropped(uint64_t xdrop,
                        std::shared_ptr<AllTips> allTips,
                        std::vector<MetagraphInterface::NodeAnnotation> & annosToBeDropped) const;
    //! remove not-well matching annos and update score of last few extension steps
    void xDrop(uint64_t xdrop,
               std::vector<std::shared_ptr<AllTips>> & tipsHis,
               bool upStream);
    //! returns the sum of scores of all annotations of all tips in current AllTips
    int totalScore() const{
        return(tipsHistory.back()->totalScore + upStreamTipsHistory.back()->totalScore);
    }
    //! all extensions steps made downstream
    std::vector<std::shared_ptr<AllTips>> tipsHistory;
    //! all extensions steps made upstream
    std::vector<std::shared_ptr<AllTips>> upStreamTipsHistory;
// private
    std::shared_ptr<MetagraphInterface const> graph;
    std::shared_ptr<Configuration const> config;
    std::shared_ptr<IdentifierMapping const> idMap;
    size_t binsize;

    // for extension analysis
    int tooManyDeletedAnnos = 0;
    int nMerges = 0;
    int nSplits = 0;
    int tooManyAnnosInInit = 0;
    int maxSteps = 0;
    int maxUpstreamSteps = 0;
};

#endif //_SeedExtension_HPP_
