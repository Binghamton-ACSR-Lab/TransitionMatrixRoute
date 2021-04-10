//
// Created by acsr on 4/5/21.
//

#ifndef TRANSITIONMATRIXROUTE_NODE_HPP
#define TRANSITIONMATRIXROUTE_NODE_HPP
#include <memory>
#include <unordered_map>
#include <eigen3/Eigen/Eigen>
#include <numeric>


namespace acsr {
    using NanowirePositionType = std::vector<std::pair<int, int>>;
    using IndexType = int64_t;
    using ControlType = uint32_t;

    NanowirePositionType indexToElectrodeVector(int n_wires, IndexType index) {
        NanowirePositionType electrode_vector(n_wires);
        for (auto i = 0; i < n_wires; ++i) {
            auto v = index & 0b1111;
            electrode_vector[n_wires - i - 1] = std::make_pair(int(v >> 2), int(v & 0b11));
            index = index >> 4;
        }
        return electrode_vector;
    }

    std::vector<IndexType> indexToSubIndexVec(IndexType index,const std::vector<int>& divided_vec) {
        std::vector<IndexType> sub_node_index(divided_vec.size());
        for(auto i=divided_vec.size()-1;i>0;--i){
            IndexType v = (1<<4*divided_vec[i])-1;
            sub_node_index[i]=index & v;
            index = index>>4*divided_vec[i];
        }
        sub_node_index[0]=(index);
        return sub_node_index;
    }

    IndexType electrodeVectorToIndex(int n_wires, const NanowirePositionType &electrode_vector) {
        IndexType index = (electrode_vector[0].first << 2) + electrode_vector[0].second;
        for (auto i = 1; i < n_wires; ++i)
            index = (index << 4) + (electrode_vector[i].first << 2) + electrode_vector[i].second;
        return index;
    }

    IndexType subIndexVecToIndex(const std::vector<IndexType>& sub_index_vec,const std::vector<int>& divided_vec) {
        auto index = sub_index_vec[0];
        for(auto i=1;i<divided_vec.size();++i)
            index = index<<(4*divided_vec[i]) | sub_index_vec[i];
        return index;
    }

    struct TransitionControlType {
        uint32_t control{};

        TransitionControlType() : control(0) {}

        explicit TransitionControlType(uint32_t _control) : control(_control) {
        }

        TransitionControlType &operator=(const TransitionControlType &myControlType) {
            control = myControlType.control;
            return *this;
        }

        TransitionControlType &operator=(uint32_t other) {
            control = other;
            return *this;
        }

        inline TransitionControlType operator+(const TransitionControlType &other) {
            return *this;
        }

        friend TransitionControlType operator+ (const TransitionControlType& lhs, const TransitionControlType &rhs) {
            return rhs;
        }
        friend bool operator== (const TransitionControlType& lhs, const TransitionControlType &rhs) {
            return lhs.control == rhs.control;
        }

        friend bool operator== (const TransitionControlType& lhs, uint32_t rhs) {
            return lhs.control == rhs;
        }
        friend bool operator< (const TransitionControlType& lhs, const TransitionControlType &rhs) {
            return lhs.control < rhs.control;
        }
        friend bool operator< (const TransitionControlType& lhs, uint32_t rhs) {
            return lhs.control < rhs;
        }


        inline TransitionControlType operator-(const TransitionControlType &other) {
            return *this;
        }



        inline TransitionControlType operator*(const TransitionControlType &other) const {
            //return TransitionControlType(0);

            if (control == 0 || other.control == 0)return TransitionControlType(0);
            auto c = control | other.control;
            auto total_control = c;
            while (c != 0) {
                if ((c & 0b11) == 0b11)return TransitionControlType(0);
                c = c >> 2;
            }
            return TransitionControlType(total_control);
        }

        inline TransitionControlType operator*(int other) const {
            if (control == 0 || other == 0)return TransitionControlType(0);
            return TransitionControlType(control);
        }


        explicit operator uint32_t() const {
            return control;
        }

        TransitionControlType &operator+=(const TransitionControlType &other) {
            control = other.control;
            return *this;
        }

        friend std::ostream &operator<<(std::ostream &out, const TransitionControlType &val) {
            out << val.control;
            return out;
        }

    };

    struct ControlVectorType{
    public:
        ControlVectorType()=default;
        explicit ControlVectorType(IndexType size):_size(size){

        }
        void setData(const std::vector<std::pair<IndexType,ControlType>>& data){
            vec = data;
            _element_size = data.size();
        }

        void setSize(IndexType size){
            _size = size;
        }

        void setData(std::vector<std::pair<IndexType,ControlType>>&& data){
            vec = std::move(data);
            _element_size = vec.size();
        }

        std::vector<IndexType> getIndexVector(){
            std::vector<IndexType> v(vec.size());
            std::transform(vec.begin(),vec.end(),v.begin(),[](const std::pair<IndexType,ControlType>& p){
                return p.first;
            });
            return v;
        }

        friend ControlVectorType operator*(const ControlVectorType& lhs,const ControlVectorType& rhs){
            ControlVectorType productor(IndexType(lhs._size*rhs._size));
            std::vector<std::pair<IndexType,ControlType>> data;
            for(auto& p1:lhs.vec){
                for(auto p2:rhs.vec){
                    if (p1.second == 0 || p2.second == 0)continue;
                    auto c = p1.second | p2.second;
                    auto total_control = c;
                    auto flag = true;
                    while (c != 0) {
                        if ((c & 0b11) == 0b11){
                            flag = false;
                            break;
                        }
                        c = c >> 2;
                    }
                    if(flag)
                    data.emplace_back(p1.first*rhs._size+p2.first,total_control);
                }
            }
            productor.setData(std::move(data));
            return productor;
        }





        std::vector<std::pair<IndexType,ControlType>> vec{};
        IndexType _size{};
        IndexType _element_size{};


    };


    class TransitionTreeNode : public std::enable_shared_from_this<TransitionTreeNode> {
    public:
        TransitionTreeNode() = default;

        TransitionTreeNode(IndexType _state, int _level) : state(_state), level(_level) {}

        virtual ~TransitionTreeNode() = default;

        TransitionTreeNode(const TransitionTreeNode &) = delete;

        TransitionTreeNode &operator=(const TransitionTreeNode &) = delete;

        std::shared_ptr<TransitionTreeNode> addChild(IndexType child_index, int _level) {
            auto child = std::make_shared<TransitionTreeNode>(child_index, _level);
            children.push_back(child);
            return child;
        }

        void addChild(const std::shared_ptr<TransitionTreeNode> &child) {
            children.push_back(child);
        }

        IndexType getState() const {
            return state;
        }

        int getLevel() const {
            return level;
        }

        void setParent(const std::shared_ptr<TransitionTreeNode> &_parent) {
            parent = _parent;
        }

        std::shared_ptr<TransitionTreeNode> &getParent() {
            return parent;
        }

        void removeChild(std::shared_ptr<TransitionTreeNode> child) {
            children.remove(child);
        }

        std::list<std::shared_ptr<TransitionTreeNode>> &getChildren() {
            return children;
        }

        void setPathQuality(int value) {
            path_quality = value;
        }

        int getPathQuality() const {
            return path_quality;
        }


    protected:
        IndexType state{};
        int level{};
        int path_quality{};
        std::shared_ptr<TransitionTreeNode> parent;
        std::list<std::shared_ptr<TransitionTreeNode>> children;

    };

    using NodePtr = std::shared_ptr<TransitionTreeNode>;
}
/*
namespace Eigen {
    template<> struct NumTraits<ACSR::TransitionControlType>
            : NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
    {
        typedef ACSR::TransitionControlType Real;
        typedef ACSR::TransitionControlType NonInteger;
        typedef ACSR::TransitionControlType Nested;
        enum {
            IsComplex = 0,
            IsInteger = 1,
            IsSigned = 0,
            RequireInitialization = 1,
            ReadCost = 1,
            AddCost = 3,
            MulCost = 3
        };
    };
    template<typename BinaryOp>
    struct ScalarBinaryOpTraits<ACSR::TransitionControlType,uint32_t,BinaryOp> { typedef ACSR::TransitionControlType ReturnType;  };

    template<typename BinaryOp>
    struct ScalarBinaryOpTraits<uint32_t ,ACSR::TransitionControlType,BinaryOp> { typedef ACSR::TransitionControlType ReturnType;  };
}

ACSR::TransitionControlType operator+(const ACSR::TransitionControlType& left,const ACSR::TransitionControlType& right){
    return right;
}*/


#endif //TRANSITIONMATRIXROUTE_NODE_HPP
