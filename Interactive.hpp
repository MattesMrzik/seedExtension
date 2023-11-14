#ifndef _Interactive_HPP_
#define _Interactive_HPP_

#include "MetagraphInterface.h"
#include "seedExtension/VisualizeGraph.hpp"

#include <iostream>

class InteractiveAndRender {
public:
    InteractiveAndRender(std::shared_ptr<Configuration const> config) {
        auto graph = config->metagraphInterface();
        std::cout << "k in Graph = " << graph->getK() << '\n';
        std::cout << "nNodes = " << graph->numNodes() << '\n';
        std::cout << "enter kmer (all capital letters) or id=int " << graph->numNodes() << " or \"show\" to render sub graph" << '\n';

        std::string s = "";
        std::cin >> s;
        size_t id = graph->getK() - 10;
        while (s != "q") {
            // read Kmer
            if(s.size() == graph->getK()) {
                id = graph->getNode(s);
            }
            // render graph
            else if((int)(s.find("show")) != -1 || (int)(s.find("Show")) != -1) {
                std::vector<MetagraphInterface::NodeID> ids{id};
                std::cout << "if current node should be rendered: enter int=pathlength" << '\n';
                std::cout << "or kmers one after the other (followed by int=pathlength) to plot multiple nodes" << '\n';
                std::string userinput = "";
                std::cin >> userinput;
                // read kmers
                //clear previous node ids
                if (userinput.size()>2) {
                    ids.clear();
                }
                // read kmers
                while (userinput.size()>2) {
                    ids.push_back(graph->getNode(userinput));
                    std::cin >> userinput;
                }
                unsigned pathlength = 1;
                try {
                    pathlength = std::stoi(userinput);
                    std::cout << "you entered: " << pathlength << '\n';
                }
                catch (std::invalid_argument & e) {
                    std::cout << "this aint a number" << '\n';
                }
                if ((int)(s.find("show")) != -1) VisualizeGraph{graph, ids, pathlength, config->binsize(), false};
                // also draw edges that dont have an extending annotation
                if ((int)(s.find("Show")) != -1) VisualizeGraph{graph, ids, pathlength, config->binsize(), true};

                int dot = system("dot -Tpng graph.gv -o graph.png");
                if (dot != 0) {
                    std::cout << "error when calling dot" << '\n';
                }
                int eog = system("eog graph.png");
                if (eog != 0) {
                    std::cout << "error when calling eog" << '\n';
                }
            }
            // node id
            else {
                try {
                    id = std::stoull(s);
                    id = id >= graph->numNodes() ? graph->numNodes() - 1 : id;
                }
                catch (std::invalid_argument & e) {
                    std::cout << "this aint a number" << '\n';
                }
            }

            std::cout << "current id = " << id << ", kmer = " << graph->getKmer(id) << '\n';
            PathBundleTip tt;
            for (auto anno : graph->getAnnotation(id)) {
                std::cout << "\t";
                tt.printMetaAnno(anno);
            }
            std::cout << '\n';
            for (auto outId : graph->getOutgoing(id)) {
                std::cout << "outgoing id = " << outId << ", kmer = " << graph->getKmer(outId) << '\n';
                for (auto anno: graph->getAnnotation(outId)) {
                    std::cout << "\t";
                    tt.printMetaAnno(anno);
                }
            }
            std::cout << '\n';
            for (auto outId : graph->getIncoming(id)) {
                std::cout << "incoming id = " << outId << ", kmer = " << graph->getKmer(outId) << '\n';
                for (auto anno: graph->getAnnotation(outId)) {
                    std::cout << "\t";
                    tt.printMetaAnno(anno);
                }
            }
            std::cin>>s;
        }
        exit(0);
    }
};

#endif // _Interactive_HPP_
