//
// Created by acsr on 4/5/21.
//

#ifndef TRANSITIONMATRIXROUTE_GLOBAL_ROUTE_HPP
#define TRANSITIONMATRIXROUTE_GLOBAL_ROUTE_HPP

#include <unordered_set>
#include <thread>
#include "node.hpp"
#include "utility.hpp"
#include "transition_matrix_constructor.hpp"
#include "boost/filesystem.hpp"
#include <mutex>
#include <SQLiteCpp/SQLiteCpp.h>

namespace acsr {
    using ControlMatrixType = std::vector<ControlVectorType>;

    class GlobalRoute {

    public:
        //constructors and deconstructor
        GlobalRoute()=default;
        virtual ~GlobalRoute();

        GlobalRoute(const GlobalRoute &) = delete;

        GlobalRoute &operator=(const GlobalRoute &) = delete;


        void init();
        void test();

        bool constructTree(int n_wire,IndexType init_index,IndexType target_index,const std::vector<int>& divided_vec);
        bool constructTree2(int n_wire, IndexType init_index, IndexType target_index, const std::vector<int> &divided_vec);
        bool constructTree3(int n_wire, IndexType init_index, IndexType target_index, const std::vector<int> &divided_vec);
        bool constructTree4(int n_wire, IndexType init_index, IndexType target_index, const std::vector<int> &divided_vec);

        std::vector<IndexType> getBestSolution();

        std::unordered_set<IndexType> transition_set;
        std::mutex transition_set_mutex;

    private:
        std::vector<ControlMatrixType> _control_matrix;
        std::vector<Eigen::SparseMatrix<int, Eigen::ColMajor, IndexType>> _transition_matrix_vec;
        std::vector<Eigen::SparseMatrix<TransitionControlType, Eigen::ColMajor, IndexType>> _control_matrix_vec;
        NodePtr root= nullptr;
        NodePtr target= nullptr;
        std::atomic_int64_t reading_times{0};

    private:
        std::vector<IndexType> getNextStepIndexVec(IndexType index,const std::vector<int>& divided_vec);
        std::vector<IndexType> getNextStepIndexVecWithinSteps(int n_wires, IndexType index, IndexType index2, const std::vector<int>& divided_vec, int step);

        std::vector<IndexType> readNextStep(IndexType index);
        void writeNextStep(IndexType,const std::vector<IndexType>& vec);

        void combineMap(std::unordered_map<IndexType,NodePtr>& m1,std::unordered_map<IndexType,NodePtr>& m2);

        template<class T1,class T2>
                Eigen::SparseVector<TransitionControlType,Eigen::ColMajor,IndexType> sparseVecMultiple(const T1& v1,const T2& v2);
        int getMinimumSteps(int n_wire,IndexType init_index,IndexType target_index,const std::vector<int>& divided_vec);
        int getMinimumSteps2(int n_wire,IndexType init_index,IndexType target_index,const std::vector<int>& divided_vec);



    };

    void GlobalRoute::init() {
        TransitionMatrixConstructor constructor;
        using namespace boost::filesystem;
        if(!exists("data") || !exists("data/transition_matrix")){
            create_directories("data/transition_matrix");
        }
        if(!exists("data/transition_matrix/_matrix_1_.bin")){
            constructor.constructMatrix(1,"data/transition_matrix/_matrix_1_.bin");
        }
        if(!exists("data/transition_matrix/_matrix_2_.bin")){
            constructor.constructMatrix(2,"data/transition_matrix/_matrix_2_.bin");
        }
        if(!exists("data/transition_matrix/_matrix_3_.bin")){
            constructor.constructMatrix(3,"data/transition_matrix/_matrix_3_.bin");
        }
        if(!exists("data/transition_matrix/_matrix_4_.bin")){
            constructor.constructMatrix(4,"data/transition_matrix/_matrix_4_.bin");
        }

        if(!exists("data/transition_matrix/_control_matrix_1_.bin")){
            constructor.constructControlMatrix(1,"data/transition_matrix/_control_matrix_1_.bin");
        }
        if(!exists("data/transition_matrix/_control_matrix_2_.bin")){
            constructor.constructControlMatrix(2,"data/transition_matrix/_control_matrix_2_.bin");
        }
        if(!exists("data/transition_matrix/_control_matrix_3_.bin")){
            constructor.constructControlMatrix(3,"data/transition_matrix/_control_matrix_3_.bin");
        }
        if(!exists("data/transition_matrix/_control_matrix_4_.bin")){
            constructor.constructControlMatrix(4,"data/transition_matrix/_control_matrix_4_.bin");
        }
        _transition_matrix_vec.resize(4);
        _control_matrix_vec.resize(4);
        _control_matrix.resize(4);
        for(auto i=0;i<4;++i){
            char ss1[200],ss2[200];
            sprintf(ss1,"data/transition_matrix/_matrix_%d_.bin",i+1);
            sprintf(ss2,"data/transition_matrix/_control_matrix_%d_.bin",i+1);
            TransitionMatrixConstructor::readTransitionMatrix(ss1,_transition_matrix_vec[i]);
            TransitionMatrixConstructor::readControlMatrix(ss2,_control_matrix_vec[i]);
            _control_matrix[i].resize(_control_matrix_vec[i].cols());
            std::vector<std::vector<std::pair<IndexType,ControlType>>> data(_control_matrix_vec[i].cols());
            for (int k=0; k<_control_matrix_vec[i].outerSize(); ++k) {
                for (Eigen::SparseMatrix<TransitionControlType, Eigen::ColMajor, IndexType>::InnerIterator it(_control_matrix_vec[i], k); it; ++it) {
                    data[it.col()].push_back(std::make_pair(it.row(),it.value().control));
                }
            }
            for(auto k=0;k<_control_matrix_vec[i].cols();++k){
                _control_matrix[i][k].setSize(_control_matrix_vec[i].rows());
                _control_matrix[i][k].setData(std::move(data[k]));
            }
        }

        /*
        if(exists("temp")){
            remove_all("temp");
        }
        create_directories("temp");
        */
        std::cout<<"Initialize Complete\n";
    }

    void GlobalRoute::test() {
        using namespace std::chrono;
        //std::random_device r;
        std::default_random_engine rd;
        std::uniform_int_distribution<int> dist(0, 3);
        auto n_wires = 8;
        long t1=0,t2=0,t3=0,t4=0;

        SQLite::Database    db("mysqlite", SQLite::OPEN_READWRITE|SQLite::OPEN_CREATE);
        db.exec("DROP TABLE IF EXISTS test");
        db.exec("CREATE TABLE test (id INTEGER PRIMARY KEY, value BLOB)");

        std::vector<IndexType> sql_data1;

        for(auto k=0;k<100;++k) {
            NanowirePositionType init_state;
            NanowirePositionType target_state;
            for (auto i = 0; i < n_wires; ++i) {
                init_state.push_back(std::make_pair(dist(rd), dist(rd)));
                target_state.push_back(std::make_pair(dist(rd), dist(rd)));
            }

            auto init_index = electrodeVectorToIndex(n_wires, init_state);

            auto start = high_resolution_clock::now();
            auto vec1 = getNextStepIndexVec(init_index, {4, 4});
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);
            t1+=duration.count();



            start = high_resolution_clock::now();
            TransitionMatrixConstructor constructor;
            auto v2 = constructor.exploreHelper(init_state);
            stop = high_resolution_clock::now();
            duration = duration_cast<microseconds>(stop - start);
            t2+=duration.count();

            //std::cout << "Total size of exploreHelper: " << v2.size() << "\n";



            /*
            start = high_resolution_clock::now();
            writeNextStep(init_index, vec1);
            //std::this_thread::sleep_for(std::chrono::milliseconds (1));
            auto v3 = readNextStep(init_index);
            stop = high_resolution_clock::now();
            duration = duration_cast<microseconds>(stop - start);
            t3+=duration.count();*/

            sql_data1.push_back(init_index);

            start = high_resolution_clock::now();
            //auto vec_size = vec1.size();


            SQLite::Statement query(db, "INSERT INTO test VALUES (?, ?)");
            query.bind(1,init_index);
            query.bind(2,vec1.data(),vec1.size()*sizeof(int64_t));
            query.exec ();
            stop = high_resolution_clock::now();
            duration = duration_cast<microseconds>(stop - start);
            t4+=duration.count();



            //std::cout << "Total size of Read From file: " << v3.size() << "\n";
        }



/*
        auto start = high_resolution_clock::now();
        for(auto v:sql_data1) {
            SQLite::Statement query(db, "SELECT * FROM test WHERE id="+std::to_string(v));

            if (query.executeStep())
            {

                SQLite::Column colSize = query.getColumn(1);
                auto vec_size = colSize.getInt();
                std::vector<int64_t> vec2(vec_size);
                SQLite::Column colBlob = query.getColumn(2);
                auto blob2 = colBlob.getBlob();
                std::memcpy(vec2.data(),blob2,vec_size*sizeof(int64_t));
            }
        }
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        t4+= duration.count();
*/
        std::cout << "Total Computation Time of getNextStepIndexVec: " << t1 << "\n";
        std::cout << "Total Computation Time of exploreHelper: " << t2 << "\n";
        std::cout << "Total Computation Time of Read From file: " << t3 << "\n";
        std::cout << "Total Computation Time of Read From Sql: " << t4 << "\n";
        std::cout<<'\n';

    }

    std::vector<IndexType> GlobalRoute::getNextStepIndexVec(IndexType index, const std::vector<int> &divided_vec) {
        /*{
            std::unique_lock<std::mutex> lock(transition_set_mutex);
            if (transition_set.find(index) != transition_set.end()) {
                ++reading_times;
                return readNextStep(index);
            }
        }*/
        auto index_vec = indexToSubIndexVec(index,divided_vec);
        std::vector<IndexType> v;
        auto size = divided_vec.size();
        if(size==1){
            return _control_matrix[divided_vec[0]-1][index_vec[0]].getIndexVector();
        }else if(size==2){
            auto return_vec = controlVectorProduct(_control_matrix[divided_vec[0]-1][index_vec[0]],_control_matrix[divided_vec[1]-1][index_vec[1]]);
            //writeNextStep(index,return_vec);
            std::unique_lock<std::mutex> lock(transition_set_mutex);
            transition_set.insert(index);
            return return_vec;
        }else{
            auto vec = _control_matrix[divided_vec[0]-1][index_vec[0]];
            for(auto i=0;i<size-1;++i){
                vec = vec*_control_matrix[divided_vec[i]-1][index_vec[i]];
            }
            auto return_vec = controlVectorProduct(vec,_control_matrix[divided_vec[size-1]-1][index_vec[size-1]]);
            //writeNextStep(index,return_vec);
            return return_vec;
        }
    }

    template<class T1, class T2>
    Eigen::SparseVector<TransitionControlType, Eigen::ColMajor, IndexType>
    GlobalRoute::sparseVecMultiple(const T1 &v1, const T2 &v2) {
        Eigen::SparseVector<TransitionControlType, Eigen::ColMajor, IndexType> vec;
        vec.resize(v1.rows()*v2.rows());
        Eigen::SparseMatrix<TransitionControlType, Eigen::ColMajor, IndexType> v = v1*v2.transpose();
        std::cout<<v1.size()<<'\n';
        std::cout<<v1.rows()<<'\n';
        std::cout<<v1.cols()<<'\n';
        std::cout<<v1.outerSize()<<'\n';
        for (int k=0; k<v.outerSize(); ++k) {
            for (Eigen::SparseMatrix<TransitionControlType, Eigen::ColMajor, IndexType>::InnerIterator it(v, k); it; ++it) {
                if(it.value()==0){
                    continue;
                }
                vec.insert(it.row()*v.cols()+it.col()) = it.value();
            }
        }
        return vec;

    }

    int GlobalRoute::getMinimumSteps2(int n_wires, IndexType init_index, IndexType target_index,
                                     const std::vector<int> &divided_vec) {


        IndexType dimension = 1;
        const int over_step = 1;
        int step = 0;
        for(auto i=0;i<n_wires;++i)
            dimension = dimension * 16;
        if(n_wires<=4){
            Eigen::SparseVector<int, Eigen::ColMajor, IndexType> sparse_state;
            sparse_state.resize(dimension);
            sparse_state.insert(init_index) = 1;
            while(sparse_state.coeffRef(target_index)==0){
                sparse_state = _transition_matrix_vec[n_wires-1]*sparse_state;
                ++step;
            }
            return step;
        }
        std::unordered_set<IndexType> index_set{init_index};
        auto min_cost = getHeuristic(n_wires,init_index,target_index);
        auto processor = std::thread::hardware_concurrency();
        std::atomic_bool final_step{false};
        std::cout<<"Try to find Minimum Steps\n";
        while(true){
            ++step;
            auto total_size = index_set.size();
            auto thread_count = std::min(IndexType (processor),IndexType(total_size));
            auto count = total_size/thread_count;
            std::vector<std::thread> thread_vec;
            std::vector<std::unordered_set<IndexType>> new_set(thread_count);
            std::cout<<"Check Step "<<step<<"; Total Size: "<<total_size<<"\n";
            for(auto i=0;i<thread_count-1;++i){
                auto& s = new_set[i];
                auto start_index = i*count;
                auto end_index = (i+1)*count;
                thread_vec.emplace_back([start_index,end_index,&index_set,this,&divided_vec,&target_index,&final_step,&s,&n_wires,&step,&min_cost](){
                    auto it = index_set.begin();
                    std::advance(it,start_index);
                    for(auto i=0;i<end_index-start_index;++i){
                        if(final_step)return;
                        auto vec = getNextStepIndexVecWithinSteps(n_wires,*it,target_index,divided_vec,min_cost + over_step-step);
                        if(step>=min_cost) {
                            for (auto v:vec) {
                                if (v == target_index){
                                    final_step = true;
                                    break;
                                }
                                s.insert(v);
                            }
                        }else{
                            for (auto v:vec) {
                                s.insert(v);
                            }
                        }
                        std::advance(it,1);
                    }
                });
            }

            auto start_index = (thread_count-1)*count;
            auto end_index = index_set.size();

            auto it = index_set.begin();
            std::advance(it,start_index);
            for(auto i=0;i<end_index-start_index;++i) {
                if (final_step)break;
                auto vec = getNextStepIndexVecWithinSteps(n_wires,*it,target_index,divided_vec,min_cost + over_step-step);
                if(step>=min_cost) {
                    for (auto v:vec) {
                        if (v == target_index){
                            final_step = true;
                            break;
                        }
                        new_set[thread_count - 1].insert(v);
                    }
                }else{
                    for (auto v:vec) {
                        new_set[thread_count - 1].insert(v);
                    }
                }
                std::advance(it, 1);
            }

            for(auto&t:thread_vec){
                if(t.joinable())t.join();
            }
            if(final_step)return step;
            index_set = std::move(new_set[0]);
            for(auto i=1;i<thread_count;++i){
                index_set.merge(std::move(new_set[i]));
            }
        }
    }

    int GlobalRoute::getMinimumSteps(int n_wires,IndexType init_index,IndexType target_index,const std::vector<int>& divided_vec) {

        IndexType dimension = 1;
        int step = 0;
        for(auto i=0;i<n_wires;++i)
            dimension = dimension * 16;
        if(n_wires<=4){
            Eigen::SparseVector<int, Eigen::ColMajor, IndexType> sparse_state;
            sparse_state.resize(dimension);
            sparse_state.insert(init_index) = 1;
            while(sparse_state.coeffRef(target_index)==0){
                sparse_state = _transition_matrix_vec[n_wires-1]*sparse_state;
                ++step;
            }
            return step;
        }
        auto min_cost = getHeuristic(n_wires,init_index,target_index);
        std::unordered_set<IndexType> index_set={init_index};
        while(true){
            ++step;
            std::unordered_set<IndexType> new_set;
            for(auto index:index_set){
                auto vec = getNextStepIndexVec(index,divided_vec);
                for(auto v:vec){
                    if(v==target_index)return step;
                    if(getHeuristic(n_wires,v,target_index)+step<min_cost+2)
                        new_set.insert(v);
                }
            }
            index_set = std::move(new_set);
        }
    }

    bool GlobalRoute::constructTree(int n_wire,IndexType init_index,IndexType target_index,const std::vector<int>& divided_vec) {
        if(init_index==target_index){
            std::cout<<"Start and Target States Are the Same!\n";
            return false;
        }
        IndexType dimension = 1;
        for(auto i=0;i<n_wire;++i)
            dimension = dimension * 16;
        target = nullptr;
        root = std::make_shared<TransitionTreeNode>(init_index,0);
        root->setPathQuality(0);

        auto start = std::chrono::high_resolution_clock::now();
        auto step = getMinimumSteps2(n_wire,init_index,target_index,divided_vec);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout<<"Find min step cost: "<<duration.count()<<"\n";
        std::vector parent_vec{root};
        for(auto i=0;i<step;++i){
            std::cout<<"Construct Tree, Step "<<i+1<<"\n";
            std::vector<std::pair<NodePtr ,std::vector<IndexType>>> children_vec;
            for(auto& parent:parent_vec){
                children_vec.emplace_back(parent,std::move(getNextStepIndexVec(parent->getState(),divided_vec)));
            }
            std::unordered_map<IndexType,NodePtr> children_map;
            for(auto& vec:children_vec){
                for(auto& child_index:vec.second){
                    if(!withinHeuristic(n_wire,target_index,child_index,step-i-1))continue;
                    auto& current_parent = children_map[child_index];
                    if( current_parent == nullptr || current_parent->getPathQuality()>vec.first->getPathQuality()){
                        current_parent = vec.first;
                    }
                }
            }
            std::vector<NodePtr> temp_parent_vec;
            for(auto& p:children_map){
                auto node = std::make_shared<TransitionTreeNode>(p.first,i+1);
                node->setPathQuality(p.second->getPathQuality()+getQuality(n_wire,p.first));
                node->setParent(p.second);
                p.second->addChild(node);
                temp_parent_vec.push_back(node);
                if(i==step-1){
                    target=node;
                }
            }
            for(auto& p:parent_vec){
                //trimLeaf(p);
            }
            parent_vec = std::move(temp_parent_vec);
        }
        if(target==nullptr){
            std::thread t(&destroyBranch,root);
            t.detach();
            return false;
        }else {
            return true;
            std::cout<<"Construct Forward Tree Successfully!\n";
        }
    }



    bool GlobalRoute::constructTree2(int n_wire,IndexType init_index,IndexType target_index,const std::vector<int>& divided_vec) {
        //transition_map.clear();
        int over_step = 1;
        if(init_index==target_index){
            std::cout<<"Start and Target States Are the Same!\n";
            return false;
        }
        IndexType dimension = 1;
        for(auto i=0;i<n_wire;++i)
            dimension = dimension * 16;


        //auto start = std::chrono::high_resolution_clock::now();
        //auto step = getMinimumSteps2(n_wire,init_index,target_index,divided_vec);

        //auto stop = std::chrono::high_resolution_clock::now();
        //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        //auto estimated_step =getHeuristic(n_wire,init_index,target_index);
        //std::cout<<"Estimated: "<<estimated_step<<"\n";
        //std::cout<<"Actual steps:"<<step<<"; Total cost: "<<duration.count()<<"\n";

        auto processor = std::thread::hardware_concurrency();
        for(over_step=0;over_step<2;++over_step) {
            auto step = getHeuristic(n_wire,init_index,target_index)+over_step;
            target = nullptr;
            root = std::make_shared<TransitionTreeNode>(init_index, 0);
            root->setPathQuality(0);
            std::vector parent_vec{root};
            for (auto i = 0; i < step; ++i) {
                auto total_size = parent_vec.size();
                auto thread_count = std::min(int(processor), int(total_size));
                auto count = total_size / thread_count;
                std::vector<std::unordered_map<IndexType, NodePtr>> children_map_vecs(thread_count);
                std::vector<std::thread> thread_vec;
                std::cout << "Construct Tree, Step " << i + 1 << "\n";
                auto start = std::chrono::high_resolution_clock::now();
                for (auto k = 0; k < thread_count - 1; ++k) {
                    auto &m = children_map_vecs[k];
                    thread_vec.emplace_back(
                            [target_index, n_wire, i, step, &parent_vec, count, k, &m, this, &divided_vec]() {
                                for (auto index = k * count; index < (k + 1) * count; ++index) {
                                    auto current_quality = parent_vec[index]->getPathQuality();
                                    auto current_parent = parent_vec[index];
                                    auto vec = getNextStepIndexVecWithinSteps(n_wire, current_parent->getState(),
                                                                              target_index, divided_vec, step - i - 1);
                                    for (auto &v:vec) {
                                        auto &p = m[v];
                                        if (p == nullptr || p->getPathQuality() > current_quality) {
                                            p = current_parent;
                                        }
                                    }
                                }
                            });
                }

                auto &m = children_map_vecs[thread_count - 1];
                for (auto index = count * (thread_count - 1); index < total_size; ++index) {
                    auto current_quality = parent_vec[index]->getPathQuality();
                    auto current_parent = parent_vec[index];
                    std::vector<IndexType> vec = getNextStepIndexVecWithinSteps(n_wire, current_parent->getState(),
                                                                                target_index, divided_vec,
                                                                                step - i - 1);
                    for (auto &v:vec) {
                        auto &p = m[v];
                        if (p == nullptr || p->getPathQuality() > current_quality) {
                            p = current_parent;
                        }
                    }
                }

                for (auto &t:thread_vec) {
                    if (t.joinable())t.join();
                }

                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
                std::cout << "Find Layer Cost: " << duration.count() << '\n';
                start = stop;


                while (thread_count > 1) {
                    auto remain = thread_count % 2;
                    int new_thread_count = thread_count / 2;
                    std::vector<std::thread> combine_thread_vec;
                    for (auto i = 0; i < new_thread_count - 1; ++i) {
                        combine_thread_vec.emplace_back(&GlobalRoute::combineMap, this, std::ref(children_map_vecs[i]),
                                                        std::ref(children_map_vecs[thread_count - 1 - i]));
                    }
                    combineMap(children_map_vecs[new_thread_count - 1],
                               children_map_vecs[thread_count - new_thread_count]);
                    for (auto &t:combine_thread_vec) {
                        if (t.joinable())t.join();
                    }
                    thread_count = new_thread_count + remain;
                }
                auto children_map = std::move(children_map_vecs[0]);


                /*
                std::unordered_map<IndexType, NodePtr> children_map;
                for (auto &children_vec:children_vecs)
                    for (auto &vec:children_vec) {
                        for (auto &child_index:vec.second) {
                            auto &current_parent = children_map[child_index];
                            if (current_parent == nullptr ||
                                current_parent->getPathQuality() > vec.first->getPathQuality()) {
                                current_parent = vec.first;
                            }
                        }
                    }*/
                stop = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
                std::cout << "Combine Cost: " << duration.count() << '\n';
                start = stop;

                std::vector<NodePtr> temp_parent_vec;
                for (auto &p:children_map) {
                    auto node = std::make_shared<TransitionTreeNode>(p.first, i + 1);
                    node->setPathQuality(p.second->getPathQuality() + getQuality(n_wire, p.first));
                    node->setParent(p.second);
                    //p.second->addChild(node);
                    temp_parent_vec.push_back(node);
                }
                if (i == step - 1 &&!temp_parent_vec.empty()) {
                    target = temp_parent_vec[0];
                }
                stop = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
                std::cout << "Add Node Cost: " << duration.count() << '\n';
                /*
                for (auto &p:parent_vec) {
                    trimLeaf(p);
                }
                 */
                parent_vec = std::move(temp_parent_vec);
            }
            if (target == nullptr) {
                std::thread t(&destroyBranch, root);
                t.detach();
            } else {
                std::cout << "Read File Times During Finding Minimum Steps Process: " << reading_times << '\n';
                std::cout << "Construct Forward Tree Successfully!\n";
                return true;
            }
        }

    }

    bool GlobalRoute::constructTree3(int n_wire,IndexType init_index,IndexType target_index,const std::vector<int>& divided_vec) {

        if(init_index==target_index){
            std::cout<<"Start and Target States Are the Same!\n";
            return false;
        }
        IndexType dimension = 1;
        for(auto i=0;i<n_wire;++i)
            dimension = dimension * 16;

        auto step = getHeuristic(n_wire,init_index,target_index)+1;
        auto processor = std::thread::hardware_concurrency();

        target= nullptr;
        root = std::make_shared<TransitionTreeNode>(init_index,0);
        root->setPathQuality(0);
        std::vector parent_vec{root};
        for (auto i = 0; i < step; ++i) {
            auto total_size = parent_vec.size();
            auto thread_count = std::min(IndexType(processor), IndexType(total_size));
            auto count = total_size / thread_count;
            std::vector<std::unordered_map<IndexType , NodePtr>> children_map_vecs(thread_count);
            std::vector<std::thread> thread_vec;
            std::cout<<"Construct Tree, Step "<<i+1<<"\n";
            auto start = std::chrono::high_resolution_clock::now();
            for (auto k = 0; k < thread_count - 1; ++k) {
                auto& m = children_map_vecs[k];
                thread_vec.emplace_back([target_index,n_wire,i,step, &parent_vec, count, k, &m, this, &divided_vec]() {
                    std::vector<std::pair<std::pair<NodePtr,int>,std::vector<IndexType>>> children_list;
                    for (auto index = k * count; index < (k + 1) * count; ++index) {
                        //auto current_quality = parent_vec[index]->getPathQuality();
                        auto current_parent = parent_vec[index];
                        auto vec = getNextStepIndexVecWithinSteps(n_wire,current_parent->getState(),target_index,divided_vec,step - i - 1);
                        children_list.push_back({{parent_vec[index],index-k * count}, std::move(vec)});
                    }
                    std::map<IndexType,std::vector<std::pair<NodePtr,int>>,std::greater<>> map;

                    for (auto it = children_list.begin(); it != children_list.end();++it) {
                        map[it->second.back()].push_back(it->first);
                        it->second.pop_back();
                    }

                    while(!map.empty()) {
                        auto& node_vec = map.begin()->second;
                        auto best_parent = node_vec.begin()->first;

                        for(auto & it : node_vec){
                            if(it.first->getPathQuality()<best_parent->getPathQuality())best_parent=it.first;
                            auto index = it.second;
                            if(!children_list[index].second.empty()){
                                map[children_list[index].second.back()].push_back(children_list[index].first);
                                children_list[index].second.pop_back();
                            }
                        }
                        m.insert(std::make_pair(map.begin()->first,best_parent));
                        map.erase(map.begin());
                    }
                });
            }

            auto& m = children_map_vecs[thread_count - 1];
            std::vector<std::pair<std::pair<NodePtr,int>,std::vector<IndexType>>> children_list;
            for (auto index = count * (thread_count - 1); index < total_size; ++index) {
                //auto current_quality = parent_vec[index]->getPathQuality();
                auto current_parent = parent_vec[index];
                std::vector<IndexType> vec = getNextStepIndexVecWithinSteps(n_wire,current_parent->getState(),target_index,divided_vec,step - i - 1);
                children_list.push_back({{parent_vec[index],index-(thread_count - 1)* count}, std::move(vec)});
            }
            std::map<IndexType,std::vector<std::pair<NodePtr,int>>,std::greater<>> map;
            std::cout<<children_list.size()<<std::endl;
            for (auto it = children_list.begin(); it != children_list.end();++it) {
                map[it->second.back()].push_back(it->first);
                it->second.pop_back();
            }

            while(!map.empty()) {
                auto& node_vec = map.begin()->second;
                auto best_parent = node_vec.begin()->first;

                for(auto & it : node_vec){
                    if(it.first->getPathQuality()<best_parent->getPathQuality())best_parent=it.first;
                    auto index = it.second;
                    if(!children_list[index].second.empty()){
                        map[children_list[index].second.back()].push_back(children_list[index].first);
                        children_list[index].second.pop_back();
                    }
                }
                m.insert(std::make_pair(map.begin()->first,best_parent));
                map.erase(map.begin());
            }


            for (auto &t:thread_vec) {
                if (t.joinable())t.join();
            }

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout<<"Find Layer Cost: "<<duration.count()<<'\n';
            start = stop;


            while(thread_count>1){
                auto remain = thread_count%2;
                int new_thread_count = thread_count/2;
                std::vector<std::thread> combine_thread_vec;
                for(auto i=0;i<new_thread_count-1;++i){
                    combine_thread_vec.emplace_back(&GlobalRoute::combineMap,this,std::ref(children_map_vecs[i]),std::ref(children_map_vecs[thread_count-1-i]));
                }
                combineMap(children_map_vecs[new_thread_count-1],children_map_vecs[thread_count-new_thread_count]);
                for(auto& t:combine_thread_vec){
                    if(t.joinable())t.join();
                }
                thread_count = new_thread_count+remain;
            }
            auto& children_map = children_map_vecs[0];


            /*
            std::unordered_map<IndexType, NodePtr> children_map;
            for (auto &children_vec:children_vecs)
                for (auto &vec:children_vec) {
                    for (auto &child_index:vec.second) {
                        auto &current_parent = children_map[child_index];
                        if (current_parent == nullptr ||
                            current_parent->getPathQuality() > vec.first->getPathQuality()) {
                            current_parent = vec.first;
                        }
                    }
                }*/
            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout<<"Combine Cost: "<<duration.count()<<'\n';
            start = stop;

            std::vector<NodePtr> temp_parent_vec;
            for (auto &p:children_map) {
                auto node = std::make_shared<TransitionTreeNode>(p.first, i + 1);
                node->setPathQuality(p.second->getPathQuality() + getQuality(n_wire, p.first));
                node->setParent(p.second);
                p.second->addChild(node);
                temp_parent_vec.push_back(node);
                if (i == step - 1) {
                    target = node;
                }
            }
            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout<<"Add Node Cost: "<<duration.count()<<'\n';
            /*
            for (auto &p:parent_vec) {
                trimLeaf(p);
            }
             */
            parent_vec = std::move(temp_parent_vec);
        }
        if (target == nullptr) {
            std::thread t(&destroyBranch, root);
            t.detach();
        } else {
            std::cout<<"Read File Times During Finding Minimum Steps Process: "<<reading_times<<'\n';
            std::cout << "Construct Forward Tree Successfully!\n";
            return true;
        }

    }

    bool GlobalRoute::constructTree4(int n_wire,IndexType init_index,IndexType target_index,const std::vector<int>& divided_vec) {

        if(init_index==target_index){
            std::cout<<"Start and Target States Are the Same!\n";
            return false;
        }
        IndexType dimension = 1;
        for(auto i=0;i<n_wire;++i)
            dimension = dimension * 16;

        auto step = getHeuristic(n_wire,init_index,target_index)+1;
        auto processor = std::thread::hardware_concurrency();

        target= nullptr;
        root = std::make_shared<TransitionTreeNode>(init_index,0);
        root->setPathQuality(0);
        std::vector parent_vec{root};
        for (auto i = 0; i < step; ++i) {
            auto total_size = parent_vec.size();
            auto thread_count = std::min(IndexType(processor), IndexType(total_size));
            auto count = total_size / thread_count;
            std::vector<std::unordered_multimap<IndexType , NodePtr>> children_map_vecs(thread_count);
            std::vector<std::thread> thread_vec;
            std::cout<<"Construct Tree, Step "<<i+1<<"\n";
            auto start = std::chrono::high_resolution_clock::now();
            for (auto k = 0; k < thread_count - 1; ++k) {
                auto& m = children_map_vecs[k];
                thread_vec.emplace_back([target_index,n_wire,i,step, &parent_vec, count, k, &m, this, &divided_vec]() {
                    for (auto index = k * count; index < (k + 1) * count; ++index) {
                        auto current_parent = parent_vec[index];
                        auto vec = getNextStepIndexVecWithinSteps(n_wire,current_parent->getState(),target_index,divided_vec,step - i - 1);
                        for(auto&v:vec){
                            m.insert({v,current_parent});
                        }
                    }
                });
            }

            auto& m = children_map_vecs[thread_count - 1];
            for (auto index = count * (thread_count - 1); index < total_size; ++index) {
                auto current_parent = parent_vec[index];
                auto vec = getNextStepIndexVecWithinSteps(n_wire,current_parent->getState(),target_index,divided_vec,step - i - 1);
                for(auto&v:vec){
                    m.insert({v,parent_vec[index]});
                }
            }

            for (auto &t:thread_vec) {
                if (t.joinable())t.join();
            }

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout<<"Find Layer Cost: "<<duration.count()<<'\n';
            start = stop;


            while(thread_count>1){
                auto remain = thread_count%2;
                int new_thread_count = thread_count/2;


                std::vector<std::thread> combine_thread_vec;
                for(auto i=0;i<new_thread_count-1;++i){
                    combine_thread_vec.emplace_back(std::thread([i,&children_map_vecs,thread_count](){
                        children_map_vecs[i].merge(children_map_vecs[thread_count-1-i]);
                        children_map_vecs[thread_count-1-i].clear();
                    }));
                }
                children_map_vecs[new_thread_count-1].merge(children_map_vecs[thread_count-new_thread_count]);
                children_map_vecs[thread_count-new_thread_count].clear();
                for(auto& t:combine_thread_vec){
                    if(t.joinable())t.join();
                }
                thread_count = new_thread_count+remain;
            }
            auto children_map = std::move(children_map_vecs[0]);

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout<<"Combine Cost: "<<duration.count()<<'\n';
            start = stop;

            std::vector<NodePtr> temp_parent_vec;
            auto best_parent = children_map.begin()->second;
            auto current_index = children_map.begin()->first;
            auto it = children_map.begin();
            std::advance(it,1);
            for (;it!=children_map.end();++it) {
                if(it->first==current_index){
                    best_parent = best_parent->getPathQuality()<=it->second->getPathQuality()?best_parent:it->second;
                    continue;
                }

                auto node = std::make_shared<TransitionTreeNode>(current_index, i + 1);
                node->setPathQuality(best_parent->getPathQuality() + getQuality(n_wire, current_index));
                node->setParent(best_parent);
                best_parent->addChild(node);
                temp_parent_vec.push_back(node);
                if (i == step - 1) {
                    target = node;
                }
            }
            auto node = std::make_shared<TransitionTreeNode>(current_index, i + 1);
            node->setPathQuality(best_parent->getPathQuality() + getQuality(n_wire, current_index));
            node->setParent(best_parent);
            best_parent->addChild(node);
            temp_parent_vec.push_back(node);
            if (i == step - 1) {
                target = node;
            }


            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout<<"Add Node Cost: "<<duration.count()<<'\n';
            /*
            for (auto &p:parent_vec) {
                trimLeaf(p);
            }
             */
            parent_vec = std::move(temp_parent_vec);
        }
        if (target == nullptr) {
            std::thread t(&destroyBranch, root);
            t.detach();
        } else {
            std::cout<<"Read File Times During Finding Minimum Steps Process: "<<reading_times<<'\n';
            std::cout << "Construct Forward Tree Successfully!\n";
            return true;
        }

    }

    GlobalRoute::~GlobalRoute() {
        boost::filesystem::remove_all("temp");
        std::thread t(&destroyBranch,root);
        t.detach();
    }

    std::vector<IndexType> GlobalRoute::getBestSolution() {
        std::vector<IndexType> solution;
        if(root== nullptr || target== nullptr)return solution;
        auto n = target;
        solution.push_back(n->getState());
        while(n->getLevel()>0){
            n = n->getParent();
            solution.push_back(n->getState());
        }
        return solution;
    }

    void GlobalRoute::writeNextStep(IndexType index, const std::vector<IndexType> &vec) {
        std::string filename = "temp/_"+std::to_string(index)+".bin";
        std::ofstream fout(filename, std::ios::out | std::ios::binary);
        fout.write((char*)&vec[0], vec.size() * sizeof(IndexType));
        fout.close();
    }

    std::vector<IndexType> GlobalRoute::readNextStep(IndexType index) {
        std::ifstream is;
        std::vector<IndexType> vec;
        std::string filename = "temp/_"+std::to_string(index)+".bin";
        is.open(filename, std::ios::binary);
        is.seekg(0, std::ios::end);
        size_t filesize=is.tellg();
        is.seekg(0, std::ios::beg);
        vec.resize(filesize/sizeof(IndexType));
        is.read((char *)vec.data(), filesize);
        return vec;
    }

    std::vector<IndexType>
    GlobalRoute::getNextStepIndexVecWithinSteps(int n_wires, IndexType index, IndexType index2, const std::vector<int> &divided_vec, int winthin_step) {
        auto index_vec = indexToSubIndexVec(index,divided_vec);
        auto size = divided_vec.size();
        if(size==1){
            std::vector<IndexType> v;
            auto vec = _control_matrix[divided_vec[0]-1][index_vec[0]].getIndexVector();
            std::copy_if(vec.begin(),vec.end(),std::back_inserter(v),[n_wires,index2,winthin_step](IndexType i){
                return withinHeuristic(n_wires, i, index2, winthin_step);
            });
            return v;
        }else{
            auto target_index_vec = indexToSubIndexVec(index2,divided_vec);
            std::vector<ControlVectorType> control_vec(size);
            for(auto i=0;i<size;++i){
                control_vec[i]=_control_matrix[divided_vec[i]-1][index_vec[i]];
            }
            return controlVectorProductWithStep(control_vec,index_vec,target_index_vec,divided_vec,winthin_step);
        }
    }

    void GlobalRoute::combineMap(std::unordered_map<IndexType, NodePtr> &map1,std::unordered_map<IndexType, NodePtr> &map2) {
        if(map1.empty() && map2.empty())return;
        for(auto & it : map2){
            auto& exist_pair = map1[it.first];
            if(exist_pair == nullptr || exist_pair->getPathQuality()>it.second->getPathQuality()){
                exist_pair = std::move(it.second);
            }
        }
        map2.clear();
    }


}
#endif //TRANSITIONMATRIXROUTE_GLOBAL_ROUTE_HPP
