#ifndef _VISUALIZEGRAPH_HPP_
#define _VISUALIZEGRAPH_HPP_

#include <memory>
#include <fstream>
#include "MetagraphInterface.h"

class VisualizeGraph {
public:
    VisualizeGraph(std::shared_ptr<MetagraphInterface const> graph_,
                    std::vector<MetagraphInterface::NodeID> nodeIDs,
                    unsigned length,
                    size_t binsize_,
                    bool drawIllegalEdges_) : graph{graph_}, binsize{binsize_}, drawIllegalEdges{drawIllegalEdges_}{
        createGVFile(nodeIDs,length);
    }
private:
    void createGVFile(std::vector<MetagraphInterface::NodeID> nodeIDs, unsigned length);

    void drawRecursive(std::ofstream & outf,
                       MetagraphInterface::NodeID nodeID,
                       unsigned length,
                       std::unordered_set<std::string>& edges,
                       bool drawIllegalEdges) const;

    void drawRecursiveInDirection(std::ofstream & outf,
                                  MetagraphInterface::NodeID nodeID,
                                  unsigned length,
                                  std::unordered_set<std::string>& edges,
                                  bool upStream,
                                  bool drawIllegalEdges) const;

    std::shared_ptr<MetagraphInterface const> graph;
    size_t binsize;
    bool drawIllegalEdges;
};


#endif //_VISUALIZEGRAPH_HPP_
