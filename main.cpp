#include <iostream>
#include <random>
#include "node.hpp"

#include "global_route.hpp"

using namespace acsr;
using namespace std::chrono;
int main() {
    auto n_wires = 8;
    std::vector<int> divided_vec{2,2,2,2};
    char s[200];
    sprintf(s,"global_route_running_time_%d_nanowires.txt",n_wires);
    std::ofstream out(s);
    //std::random_device rd;
    std::default_random_engine rd;
    std::uniform_int_distribution<int> dist(0, 3);

    for(auto k=0;k<20;++k) {
        GlobalRoute globalRoute;
        //globalRoute.readTransitiomMatrices();
        globalRoute.init();
        std::cout<<"Iteration: "<<k<<std::endl;
        NanowirePositionType init_state;
        NanowirePositionType target_state;
        for (auto i = 0; i < n_wires; ++i) {
            init_state.push_back(std::make_pair(dist(rd),dist(rd)));
            target_state.push_back(std::make_pair(dist(rd),dist(rd)));
        }
        auto init_index = electrodeVectorToIndex(n_wires,init_state);
        auto target_index = electrodeVectorToIndex(n_wires,target_state);

        auto start = high_resolution_clock::now();
        std::vector<IndexType> path;


        if(globalRoute.constructTree2(n_wires,init_index,target_index,divided_vec)){
            path = globalRoute.getBestSolution();
            for(auto i :path)std::cout<<i<<"->";
            std::cout<<"\n";
        }else{
            std::cout<<"cannot find solution\n";
            continue;
        }

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        out<<duration.count()<<"\n";

        exportSVG(n_wires,k,path);
        auto estimate_steps = getHeuristic(n_wires,init_state,target_state);
        writeToDatabase(n_wires,init_state,target_state,estimate_steps,path.size()-estimate_steps-1,divided_vec,duration.count(),path);

    }
    out.close();

    /*
    using namespace ACSR;
    using namespace std::chrono;

    //std::random_device r;
    std::default_random_engine rd;
    std::uniform_int_distribution<int> dist(0, 3);

    ACSR::GlobalRoute route;
    route.init();
    //route.test();
    auto n_wires = 8;

    NanowirePositionType init_state={{1,3},{0,3},{3,1},{0,3},{0,3},{2,3},{3,3},{2,1}  };
    NanowirePositionType target_state={{0,1},{1,2},{1,2},{3,3},{2,0},{2,2},{1,0},{3,2}  };
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
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout<<"Total Computation Time: "<<duration.count()<<"\n";

    return 0;*/
}
