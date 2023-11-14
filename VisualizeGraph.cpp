#include "VisualizeGraph.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <memory>

#include "MetagraphInterface.h"

#include "PathBundleTip.hpp"

void VisualizeGraph::createGVFile(std::vector<MetagraphInterface::NodeID> nodeIDs, unsigned length) {
    std::ofstream outf;
    outf.open("graph.gv");
    outf << "digraph {\n rankdir=LR;\n";
    std::unordered_set<std::string> edges;
    for (auto nodeID : nodeIDs) {
        drawRecursive(outf, nodeID, length, edges, drawIllegalEdges);
        //TODO split also init kmer
        // if (leftKmer.size() > 20) {
        //     unsigned middle = leftKmer.size()/2;
        //     leftKmerSplitted = leftKmer.substr(0,middle) + "\\n" + leftKmer.substr(middle + 1, leftKmer.size() - 1);
        //     rightKmerSplitted = rightKmer.substr(0,middle) + "\\n" + rightKmer.substr(middle + 1, rightKmer.size() - 1);
        // }
        outf << "\"" << graph->getKmer(nodeID) << "\\n" <<nodeID << "\"" << "[shape=rectangle];\n";
    }
    outf << "}\n";
    outf.close();
}

void VisualizeGraph::drawRecursive(std::ofstream & outf,
                                   MetagraphInterface::NodeID nodeID,
                                   unsigned length,
                                   std::unordered_set<std::string>& edges,
                                   bool drawIllegalEdges) const {
    if (length > 0) {
        drawRecursiveInDirection(outf, nodeID, length, edges, false, drawIllegalEdges);
        drawRecursiveInDirection(outf, nodeID, length, edges, true, drawIllegalEdges);
    }
}

void VisualizeGraph::drawRecursiveInDirection(std::ofstream & outf,
                                              MetagraphInterface::NodeID nodeID,
                                              unsigned length,
                                              std::unordered_set<std::string>& edges,
                                              bool upStream,
                                              bool drawIllegalEdges) const{// edges, where the intersection of annotations of both nodes is empty
    std::string currentKmer = graph->getKmer(nodeID);
    auto && currentAnnotations = graph->getAnnotation(nodeID);
    auto outIDs = upStream ? graph->getIncoming(nodeID) : graph->getOutgoing(nodeID);
    for (auto && outID : outIDs) {
        bool foundMatchingAnno = false;
        auto && outgoingKmer = graph->getKmer(outID);
        std::string label = "[label=\"";
        auto && nextAnnotations = graph->getAnnotation(outID);

        //search for matchin anno
        for (auto & metaAnno : currentAnnotations) {
            for (auto & nextMetaAnno : nextAnnotations) {
                bool isTransitionInRightDirection = upStream ?
                                                    metaAnno.bin_idx >= nextMetaAnno.bin_idx :
                                                    metaAnno.bin_idx <= nextMetaAnno.bin_idx;
                if (metaAnno.genome == nextMetaAnno.genome && //found matching genome
                    metaAnno.sequence == nextMetaAnno.sequence &&
                    isTransitionInRightDirection &&
                    labs(metaAnno.bin_idx - nextMetaAnno.bin_idx) < 2 * binsize) { //step isnt to big

                    foundMatchingAnno = true;
                    //bin_idx of anno, with possible transition
                    auto leftBinIdx  = upStream ? std::to_string(nextMetaAnno.bin_idx/binsize) : std::to_string(metaAnno.bin_idx/binsize);
                    auto rightBinIdx = upStream ? std::to_string(metaAnno.bin_idx/binsize) : std::to_string(nextMetaAnno.bin_idx/binsize);
                    auto coord = metaAnno.bin_idx == nextMetaAnno.bin_idx ?
                                 leftBinIdx :
                                 leftBinIdx + " -> " + rightBinIdx;
                    label.append(metaAnno.genome.substr(0, metaAnno.genome.size() - 3) + ", "
                                 + coord + ",\\n");
                }
            }
        }
        auto leftKmer  = upStream ? outgoingKmer : currentKmer;
        auto rightKmer = upStream ? currentKmer : outgoingKmer;
        auto leftKmerSplitted  = upStream ? outgoingKmer : currentKmer;
        auto rightKmerSplitted = upStream ? currentKmer : outgoingKmer;
        // 
        // if (leftKmer.size() > 20) {
        //     unsigned middle = leftKmer.size()/2;
        //     leftKmerSplitted = leftKmer.substr(0,middle) + "\\n" + leftKmer.substr(middle + 1, leftKmer.size() - 1);
        //     rightKmerSplitted = rightKmer.substr(0,middle) + "\\n" + rightKmer.substr(middle + 1, rightKmer.size() - 1);
        // }

        bool a = foundMatchingAnno;
        bool b = edges.find(currentKmer + outgoingKmer) == edges.end();
        bool c = drawIllegalEdges;
        //I designed a truth table with desired output and applied online alg
        //http://www.32x8.com/sop3_____A-B-C_____m_1-4-5___________option-0_____999790975275833392714
        if ((a && b) || (b && c)) {
            //draw edge with both nodes
            outf << "\""
            << leftKmerSplitted << "\\n" << graph->getNode(leftKmer)
            << "\" -> \""
            << rightKmerSplitted << "\\n" << graph->getNode(rightKmer)
            << "\"";
            //draw label
            outf << label << "\"];\n";
            edges.insert(currentKmer + outgoingKmer);
            edges.insert(outgoingKmer + currentKmer);

            drawRecursive(outf, outID, length - 1, edges, drawIllegalEdges);
        }
    }
}
