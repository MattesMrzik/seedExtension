#ifndef _PATHBUNDLETIP_HPP_
#define _PATHBUNDLETIP_HPP_

#include "MetagraphInterface.h"

#include <vector>
#include <memory>


/*! Represents all annotations of a metagraph node
* which are considered in the current seed extension.
* \details Each annotation is a tuple of metagraph annotation and Score.
* Every annotation a_ has its own score, since we want to drop a_ in xDrop
* if a_ doenst match well to all the other annotations
* (Knoten v, Buendel B(v))
*/
class PathBundleTip{

public:
    struct HashNodeAnnotation{
        std::size_t operator()(MetagraphInterface::NodeAnnotation const & anno) const {
            return (std::hash<std::string>()(anno.genome)
                    ^ std::hash<std::string>()(anno.sequence)
                    ^ std::hash<bool>()(anno.reverse_strand)
                    ^ std::hash<int>()(anno.bin_idx));
        }
    };
    //this is the value associated to a metagraph annotation
    struct Score{
        // current score of annotation
        double currentScore;
        // the max score this annotation had at some point of the seed extension process
        double maxScore;
        // in which extension step was the max score set
        uint64_t ageOfMaxScore;
        // at what extension step was the latest bin_idx transition made
        uint64_t latestTransition;
    };

    struct HashPathBundleTip{
        std::size_t operator()(uint64_t const & nodeID) const {
            return(nodeID);
        }
    };

    PathBundleTip(uint64_t nodeId_,
                  std::unordered_map<MetagraphInterface::NodeAnnotation,
                                     Score,
                                     HashNodeAnnotation> annotations_):
                  nodeID{nodeId_},
                  annotations{annotations_}{}

    PathBundleTip():nodeID{0},annotations{}{}

    /*! searches in adjacent nodes, whether the sequence corresponding to
    * a current annotation continues
    */
    std::vector<std::shared_ptr<PathBundleTip>> extendTip(std::shared_ptr<MetagraphInterface const>  graph,
                                         bool upStream,
                                         size_t binsize,
                                         uint64_t numberOfExtensionsMade);

    void print(std::shared_ptr<MetagraphInterface const>  graph) const;

    void print() const; // without kmer

    static void printMetaAnno(MetagraphInterface::NodeAnnotation const & anno);

    bool operator==(PathBundleTip const & other) const{
        return (nodeID == other.nodeID);
    }

    uint64_t nodeID;
    // all annotations for this bundle
    std::unordered_map<MetagraphInterface::NodeAnnotation,
                       Score,
                       HashNodeAnnotation> annotations;

};
#endif //_PATHBUNDLETIP_HPP_
