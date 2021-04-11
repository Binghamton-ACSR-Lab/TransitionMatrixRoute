#include <iostream>
#include <random>
#include "node.hpp"

#include "global_route.hpp"

using namespace acsr;
using namespace std::chrono;
int main() {
    auto n_wires = 8;
    std::vector<int> divided_vec{2,2,2,2};

    //std::random_device rd;
    std::default_random_engine rd;
    std::uniform_int_distribution<int> dist(0, 3);

    for(auto k=0;k<20;++k) {
        GlobalRoute globalRoute;
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


        if(globalRoute.constructTree(n_wires,init_index,target_index,divided_vec)){
            path = globalRoute.getBestSolution();
            for(auto i :path)std::cout<<i<<"->";
            std::cout<<"\n";
        }else{
            std::cout<<"cannot find solution\n";
            continue;
        }

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);


        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
        exportSVG(n_wires,path,"image/image_"+ss.str()+".svg");
        auto estimate_steps = getHeuristic(n_wires,init_state,target_state);
        writeToDatabase(n_wires,init_state,target_state,estimate_steps,path.size()-estimate_steps-1,divided_vec,duration.count(),path);

    }

}
