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

namespace ACSR {
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

        std::vector<IndexType> getBestSolution();

        std::unordered_map<IndexType,std::vector<IndexType>> transition_map;
        std::mutex transition_map_mutex;

    private:
        std::vector<ControlMatrixType> _control_matrix;
        std::vector<Eigen::SparseMatrix<int, Eigen::ColMajor, IndexType>> _transition_matrix_vec;
        std::vector<Eigen::SparseMatrix<TransitionControlType, Eigen::ColMajor, IndexType>> _control_matrix_vec;

    private:
        std::vector<IndexType> getNextStepIndexVec(IndexType index,const std::vector<int>& divided_vec);

        template<class T1,class T2>
                Eigen::SparseVector<TransitionControlType,Eigen::ColMajor,IndexType> sparseVecMultiple(const T1& v1,const T2& v2);
        int getMinimumSteps(int n_wire,IndexType init_index,IndexType target_index,const std::vector<int>& divided_vec);
        int getMinimumSteps2(int n_wire,IndexType init_index,IndexType target_index,const std::vector<int>& divided_vec);

        NodePtr root= nullptr;
        NodePtr target= nullptr;


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
            ACSR::TransitionMatrixConstructor::readTransitionMatrix(ss1,_transition_matrix_vec[i]);
            ACSR::TransitionMatrixConstructor::readControlMatrix(ss2,_control_matrix_vec[i]);
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

        std::cout<<"Initialize Complete\n";
    }

    void GlobalRoute::test() {
        using namespace std::chrono;
        std::random_device r;
        std::default_random_engine rd(r());
        std::uniform_int_distribution<int> dist(0, 3);
        auto n_wires = 8;
        NanowirePositionType init_state;
        NanowirePositionType target_state;
        for (auto i = 0; i < n_wires; ++i) {
            init_state.push_back(std::make_pair(dist(rd),dist(rd)));
            target_state.push_back(std::make_pair(dist(rd),dist(rd)));
        }

        auto init_index = electrodeVectorToIndex(n_wires,init_state);
        auto start = high_resolution_clock::now();

        auto vec1 = getNextStepIndexVec(init_index,{4,4});

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        std::cout<<"Total Computation Time of getNextStepIndexVec: "<<duration.count()<<"\n";
        std::cout<<"Total size of getNextStepIndexVec: "<<vec1.size()<<"\n";
        start = stop;
        TransitionMatrixConstructor constructor;
        auto v2 = constructor.exploreHelper(init_state);
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - start);
        std::cout<<"Total Computation Time of exploreHelper: "<<duration.count()<<"\n";
        std::cout<<"Total size of exploreHelper: "<<v2.size()<<"\n";
        std::vector<IndexType> vec2(v2.size());
        std::transform(v2.begin(),v2.end(),vec2.begin(),[n_wires](const NanowirePositionType& p){
            return electrodeVectorToIndex(n_wires,p);
        });
        std::sort(vec2.begin(),vec2.end());

        auto s = v2.size();
        for(auto i=0;i<s;++i){
            if(vec1[i]!=vec2[i]){
                std::cout<<vec1[i]<<std::endl;
            }
        }
        std::cout<<'\n';

    }

    std::vector<IndexType> GlobalRoute::getNextStepIndexVec(IndexType index, const std::vector<int> &divided_vec) {
        if(transition_map.find(index)!=transition_map.end())
            return transition_map[index];
        auto index_vec = indexToSubIndexVec(index,divided_vec);
        std::vector<IndexType> v;
        auto size = divided_vec.size();
        if(size==1){
            return _control_matrix[divided_vec[0]-1][index_vec[0]].getIndexVector();
        }else if(size==2){
            auto return_vec = controlVectorProduct(_control_matrix[divided_vec[0]-1][index_vec[0]],_control_matrix[divided_vec[1]-1][index_vec[1]]);
            std::unique_lock<std::mutex> lock(transition_map_mutex);
            transition_map[index] = return_vec;
            return return_vec;
        }else{
            auto vec = _control_matrix[divided_vec[0]-1][index_vec[0]];
            for(auto i=0;i<size-1;++i){
                vec = vec*_control_matrix[divided_vec[i]-1][index_vec[i]];
            }
            auto return_vec = controlVectorProduct(vec,_control_matrix[divided_vec[size-1]-1][index_vec[size-1]]);
            std::unique_lock<std::mutex> lock(transition_map_mutex);
            transition_map[index] = return_vec;
            return return_vec;
        }

        /*
        if(divided_vec.size()==1){
            Eigen::SparseVector<TransitionControlType, Eigen::ColMajor, IndexType> vec = _control_matrix_vec[divided_vec[0]-1].col(index_vec[0]);
            for (Eigen::SparseVector<TransitionControlType, Eigen::ColMajor, IndexType>::InnerIterator it(vec); it; ++it) {
                v.push_back(it.row());
            }
        }else{
            Eigen::SparseVector<TransitionControlType, Eigen::ColMajor, IndexType> vec = _control_matrix_vec[divided_vec[0]-1].col(index_vec[0]);
            for(auto i=1;i<divided_vec.size()-1;++i){
                vec = sparseVecMultiple(vec,_control_matrix_vec[divided_vec[i]-1].col(index_vec[i]));
            }
            Eigen::SparseMatrix<TransitionControlType, Eigen::ColMajor, IndexType> m = vec*_control_matrix_vec[divided_vec.back()-1].col(index_vec.back()).transpose();
            for (int k=0; k<m.outerSize(); ++k) {
                for (Eigen::SparseMatrix<TransitionControlType, Eigen::ColMajor, IndexType>::InnerIterator it(m, k); it; ++it) {
                    if(it.value()==0){
                        continue;
                    }
                    v.push_back(it.row()*m.cols()+it.col());
                }
            }
            return v;

        }*/
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
        const int over_step = 2;
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
        while(true){
            ++step;
            auto total_size = index_set.size();
            auto thread_count = std::min(IndexType (processor),IndexType(total_size));
            auto count = total_size/thread_count;
            std::vector<std::thread> thread_vec;
            std::vector<std::unordered_set<IndexType>> new_set(thread_count);
            for(auto i=0;i<thread_count-1;++i){
                auto& s = new_set[i];
                auto start_index = i*count;
                auto end_index = (i+1)*count;
                thread_vec.emplace_back([start_index,end_index,&index_set,this,&divided_vec,&target_index,&final_step,&s,&n_wires,&step,&min_cost](){
                    auto it = index_set.begin();
                    std::advance(it,start_index);
                    for(auto i=0;i<end_index-start_index;++i){
                        if(final_step)return;
                        auto vec = getNextStepIndexVec(*it,divided_vec);
                        for(auto v:vec){
                            if(v==target_index)final_step = true;
                            if(getHeuristic(n_wires,v,target_index)+step<min_cost+over_step)
                                s.insert(v);
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
                auto vec = getNextStepIndexVec(*it, divided_vec);
                for (auto v:vec) {
                    if (v == target_index)final_step = true;
                    if(getHeuristic(n_wires,v,target_index)+step<min_cost+over_step)
                        new_set[thread_count-1].insert(v);
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
            std::vector<std::pair<NodePtr ,std::vector<IndexType>>> children_vec;
            for(auto& parent:parent_vec){
                children_vec.emplace_back(parent,std::move(getNextStepIndexVec(parent->getState(),divided_vec)));
            }
            std::unordered_map<IndexType,NodePtr> children_map;
            for(auto& vec:children_vec){
                for(auto& child_index:vec.second){
                    if(getHeuristic(n_wire,target_index,child_index)>step-i-1)continue;
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
                trimLeaf(p);
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
        transition_map.clear();
        if(init_index==target_index){
            std::cout<<"Start and Target States Are the Same!\n";
            return false;
        }
        IndexType dimension = 1;
        for(auto i=0;i<n_wire;++i)
            dimension = dimension * 16;


        auto start = std::chrono::high_resolution_clock::now();
        auto step = getMinimumSteps2(n_wire,init_index,target_index,divided_vec);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        auto estimated_step =getHeuristic(n_wire,init_index,target_index);
        std::cout<<"Estimated: "<<estimated_step<<"\n";
        std::cout<<"Actual steps:"<<step<<"; Total cost: "<<duration.count()<<"\n";

        auto processor = std::thread::hardware_concurrency();

            target= nullptr;
            root = std::make_shared<TransitionTreeNode>(init_index,0);
            root->setPathQuality(0);
            std::vector parent_vec{root};
            for (auto i = 0; i < step; ++i) {
                auto total_size = parent_vec.size();
                auto thread_count = std::min(IndexType(processor), IndexType(total_size));
                auto count = total_size / thread_count;
                std::vector<std::vector<std::pair<NodePtr, std::vector<IndexType>>>> children_vecs(thread_count);
                std::vector<std::thread> thread_vec;
                auto f = [&n_wire, &target_index, &step, i](IndexType x) {
                    return getHeuristic(n_wire, target_index, x) <= step - i - 1;
                };
                for (auto k = 0; k < thread_count - 1; ++k) {
                    thread_vec.emplace_back([&f, &parent_vec, count, k, &children_vecs, this, &divided_vec]() {
                        for (auto index = k * count; index < (k + 1) * count; ++index) {
                            auto vec = getNextStepIndexVec(parent_vec[index]->getState(), divided_vec);
                            std::vector<IndexType> v;
                            std::copy_if(vec.begin(), vec.end(), std::back_inserter(v), [&f](IndexType x) {
                                return f(x);
                            });
                            children_vecs[k].emplace_back(parent_vec[index], std::move(v));
                        }
                    });
                }

                for (auto index = count * (thread_count - 1); index < total_size; ++index) {
                    auto vec = getNextStepIndexVec(parent_vec[index]->getState(), divided_vec);
                    std::vector<IndexType> v;
                    std::copy_if(vec.begin(), vec.end(), std::back_inserter(v), [&f](IndexType x) {
                        return f(x);
                    });
                    children_vecs[thread_count - 1].emplace_back(parent_vec[index], std::move(v));
                }

                for (auto &t:thread_vec) {
                    if (t.joinable())t.join();
                }


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
                    }
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
                for (auto &p:parent_vec) {
                    trimLeaf(p);
                }
                parent_vec = std::move(temp_parent_vec);
            }
            if (target == nullptr) {
                std::thread t(&destroyBranch, root);
                t.detach();
            } else {
                std::cout << "Construct Forward Tree Successfully!\n";
                return true;
            }

    }

    GlobalRoute::~GlobalRoute() {
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


}
#endif //TRANSITIONMATRIXROUTE_GLOBAL_ROUTE_HPP
