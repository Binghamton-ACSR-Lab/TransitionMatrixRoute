#include <iostream>
#include <random>
#include "node.hpp"

#include "global_route.hpp"


int main() {
    using namespace ACSR;
    using namespace std::chrono;

    std::random_device r;
    std::default_random_engine rd(r());
    std::uniform_int_distribution<int> dist(0, 3);

    ACSR::GlobalRoute route;
    route.init();
    route.test();
    auto n_wires = 8;
    NanowirePositionType init_state;
    NanowirePositionType target_state;
    for (auto i = 0; i < n_wires; ++i) {
        init_state.push_back(std::make_pair(dist(rd),dist(rd)));
        target_state.push_back(std::make_pair(dist(rd),dist(rd)));
    }
    auto start = high_resolution_clock::now();
    std::vector<IndexType> path;
    std::pair<NodePtr,NodePtr> tree;
    std::vector<int> divided_vec{4,4};

    auto init_index = electrodeVectorToIndex(n_wires,init_state);
    auto target_index = electrodeVectorToIndex(n_wires,target_state);
    if(route.constructTree2(n_wires,init_index,target_index,divided_vec)){
        auto path = route.getBestSolution();
        for(auto&p:path){
            std::cout<<p<<"->";
        }
        std::cout<<'\n';
    };

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout<<"Total Computation Time: "<<duration.count()<<"\n";

    return 0;
}
