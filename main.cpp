#include <iostream>
#include <random>
#include "node.hpp"

#include "global_route.hpp"
#include "acsr_astar.hpp"

using namespace acsr;
using namespace std::chrono;

std::vector<std::vector<int>> divided_vecs{{2},{3},{4},{3,2},{3,3},{4,3},{4,4},{3,3,3},{3,3,4},{3,4,4},{4,4,4},{4,3,3,3},{4,4,3,3}};

void transitionVectorCompare(int n_wires){
    std::cout<<"transitionVectorCompare\nWire Count = "<<n_wires<<std::endl;
    long transition_time = 0;
    long explorer_time = 0;
    std::default_random_engine rd;
    std::uniform_int_distribution<int> dist(0, 3);

    TransitionMatrixConstructor constructor;
    GlobalRoute globalRoute;
    globalRoute.init();

    SQLite::Database db("plannerDB.db",SQLite::OPEN_READWRITE|SQLite::OPEN_CREATE);
    {
        SQLite::Statement query(db, "create table if not exists Vector ("
                                    "id INTEGER PRIMARY KEY AUTOINCREMENT, "
                                    "nanowire_count INTEGER NOT NULL,"
                                    "transition_time INTEGER NOT NULL,"
                                    "explore_solution_time INTEGER NOT NULL"
                                    ")");
        try {
            query.exec();
        } catch (SQLite::Exception &e) {
            std::cout << "Create Table Error\n";
            std::cout << e.what();
            return;
        }
    }

    for (auto k = 0; k < 30; ++k) {
        std::cout << "Iteration: " << k << std::endl;
        NanowirePositionType init_state;
        for (auto i = 0; i < n_wires; ++i) {
            init_state.push_back(std::make_pair(dist(rd), dist(rd)));
        }
        auto init_index = electrodeVectorToIndex(n_wires, init_state);

        auto start = high_resolution_clock::now();
        globalRoute.getNextStepIndexVec(n_wires,init_index,divided_vecs[n_wires-2]);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        transition_time=duration.count();

        start = high_resolution_clock::now();
        constructor.exploreHelper(init_state);
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        explorer_time=duration.count();


        {
            std::string query_string = "INSERT INTO Vector (nanowire_count,transition_time,explore_solution_time) VALUES(?,?,?)";
            SQLite::Statement query(db, query_string);


            query.bind(1, n_wires);
            query.bind(2, transition_time);
            query.bind(3, explorer_time);

            try {
                query.exec();
            } catch (SQLite::Exception &e) {
                std::cout << "Insert Solution Error\n";
                std::cout << e.what();
                return;
            }
        }

    }



}

void runTransitionMatrixPlanner(int n_wires){
    std::default_random_engine rd;
    std::uniform_int_distribution<int> dist(0, 3);
    for (auto k = 0; k < 30; ++k) {
        GlobalRoute globalRoute;
        globalRoute.init();
        std::cout << "Iteration: " << k << std::endl;
        NanowirePositionType init_state;
        NanowirePositionType target_state;
        for (auto i = 0; i < n_wires; ++i) {
            init_state.push_back(std::make_pair(dist(rd), dist(rd)));
            target_state.push_back(std::make_pair(dist(rd), dist(rd)));
        }
        auto init_index = electrodeVectorToIndex(n_wires, init_state);
        auto target_index = electrodeVectorToIndex(n_wires, target_state);
        std::vector<IndexType> path;

        auto start = high_resolution_clock::now();

        if(globalRoute.constructTree(n_wires,init_index,target_index,divided_vecs[n_wires-2])){
            path = globalRoute.getBestSolution();
            for(auto i :path)std::cout<<i<<"->";
            std::cout<<"\n";
        }else{
            std::cout<<"cannot find solution\n";
            continue;
        }

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");

        auto estimate_steps = getHeuristic(n_wires,init_state,target_state);
        std::vector<long> qualities(path.size());
        std::transform(path.begin(),path.end(),qualities.begin(),[n_wires](IndexType index){
            return -getQuality(n_wires,index);
        });
        auto path_quality = std::accumulate(qualities.begin(),qualities.end(),0L);

        writeToDatabase(n_wires,init_state,target_state,estimate_steps,path.size()-estimate_steps-1,divided_vecs[n_wires-2],duration.count(),path_quality,path);

    }
}

void runAStarPlanner(int n_wires) {
    std::default_random_engine rd;
    std::uniform_int_distribution<int> dist(0, 3);
    for (auto k = 0; k < 30; ++k) {
        GlobalRoute globalRoute;
        globalRoute.init();
        std::cout << "Iteration: " << k << std::endl;
        NanowirePositionType init_state;
        NanowirePositionType target_state;
        for (auto i = 0; i < n_wires; ++i) {
            init_state.push_back(std::make_pair(dist(rd), dist(rd)));
            target_state.push_back(std::make_pair(dist(rd), dist(rd)));
        }
        auto init_index = electrodeVectorToIndex(n_wires, init_state);
        auto target_index = electrodeVectorToIndex(n_wires, target_state);

        //auto start = high_resolution_clock::now();


        AcsrAstar astar;
        astar.init(n_wires, init_index, target_index, divided_vecs[n_wires - 2]);
        astar.run();
        //auto stop = high_resolution_clock::now();
        //auto duration = duration_cast<milliseconds>(stop - start);

        auto best_path = astar.getBestSolution();
        auto init_path = astar.init_path;

        std::vector<long> best_qualities(best_path.size());
        std::transform(best_path.begin(), best_path.end(), best_qualities.begin(), [n_wires](IndexType index) {
            return -getQuality(n_wires, index);
        });
        auto best_path_quality = std::accumulate(best_qualities.begin(), best_qualities.end(), 0L);

        std::vector<long> init_qualities(init_path.size());
        std::transform(init_path.begin(), init_path.end(), init_qualities.begin(), [n_wires](IndexType index) {
            return -getQuality(n_wires, index);
        });
        auto init_path_quality = std::accumulate(init_qualities.begin(), init_qualities.end(), 0L);

        writeToDatabaseAStar(n_wires, init_state, target_state, divided_vecs[n_wires - 2], astar.getFirstSolutionTime(),
                             astar.getBestSolutionTime(), astar.getTotalRunningTime(), init_path.size() - 1,
                             init_path_quality, best_path.size() - 1, best_path_quality, best_path);


    }
}


int main() {
    for(auto n_wires = 2;n_wires<=8;++n_wires) {
        runTransitionMatrixPlanner(n_wires);
    }
}




