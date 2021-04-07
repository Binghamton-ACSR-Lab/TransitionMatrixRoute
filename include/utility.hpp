//
// Created by acsr on 4/5/21.
//

#ifndef TRANSITIONMATRIXROUTE_UTILITY_HPP
#define TRANSITIONMATRIXROUTE_UTILITY_HPP
#include "node.hpp"
#include <fstream>
namespace ACSR {
    enum SteerDirection {
        Up = 0,
        Down,
        Left,
        Right,
        Stay
    };

    template<class SparseMatrix>
    inline void writeSparseMatrixToBin(const std::string &filename, const SparseMatrix &matrix) {
        assert(matrix.isCompressed() == true);
        std::ofstream out(filename, std::ios::binary | std::ios::out | std::ios::trunc);
        if (out.is_open()) {
            typename SparseMatrix::Index rows, cols, nnzs, outS, innS;
            rows = matrix.rows();
            cols = matrix.cols();
            nnzs = matrix.nonZeros();
            outS = matrix.outerSize();
            innS = matrix.innerSize();

            out.write(reinterpret_cast<char *>(&rows), sizeof(typename SparseMatrix::Index));
            out.write(reinterpret_cast<char *>(&cols), sizeof(typename SparseMatrix::Index));
            out.write(reinterpret_cast<char *>(&nnzs), sizeof(typename SparseMatrix::Index));
            out.write(reinterpret_cast<char *>(&outS), sizeof(typename SparseMatrix::Index));
            out.write(reinterpret_cast<char *>(&innS), sizeof(typename SparseMatrix::Index));

            typename SparseMatrix::Index sizeIndexS = static_cast<typename SparseMatrix::Index>(sizeof(typename SparseMatrix::Index));
            typename SparseMatrix::Index sizeScalar = static_cast<typename SparseMatrix::Index>(sizeof(typename SparseMatrix::Scalar));
            out.write(reinterpret_cast<const char *>(matrix.valuePtr()), sizeScalar * nnzs);
            out.write(reinterpret_cast<const char *>(matrix.outerIndexPtr()), sizeIndexS * outS);
            out.write(reinterpret_cast<const char *>(matrix.innerIndexPtr()), sizeIndexS * nnzs);
            out.close();
        } else {
            std::cout << "Can not write to file: " << filename << std::endl;
        }
    }

    template<class SparseMatrix>
    inline void readSparsMatrixFromBin(const std::string &filename, SparseMatrix &matrix) {
        std::ifstream in(filename, std::ios::binary | std::ios::in);
        if (in.is_open()) {
            typename SparseMatrix::Index rows, cols, nnz, inSz, outSz;
            typename SparseMatrix::Index sizeScalar = static_cast<typename SparseMatrix::Index>(sizeof(typename SparseMatrix::Scalar));
            typename SparseMatrix::Index sizeIndex = static_cast<typename SparseMatrix::Index>(sizeof(typename SparseMatrix::Index));
            typename SparseMatrix::Index sizeIndexS = static_cast<typename SparseMatrix::Index>(sizeof(typename SparseMatrix::StorageIndex));
            std::cout << sizeScalar << " " << sizeIndex << std::endl;
            in.read(reinterpret_cast<char *>(&rows ), sizeIndex);
            in.read(reinterpret_cast<char *>(&cols ), sizeIndex);
            in.read(reinterpret_cast<char *>(&nnz  ), sizeIndex);
            in.read(reinterpret_cast<char *>(&outSz), sizeIndex);
            in.read(reinterpret_cast<char *>(&inSz ), sizeIndex);

            matrix.resize(rows, cols);
            matrix.makeCompressed();
            matrix.resizeNonZeros(nnz);

            in.read(reinterpret_cast<char *>(matrix.valuePtr()), sizeScalar * nnz);
            in.read(reinterpret_cast<char *>(matrix.outerIndexPtr()), sizeIndexS * outSz);
            in.read(reinterpret_cast<char *>(matrix.innerIndexPtr()), sizeIndexS * nnz);

            matrix.finalize();
            in.close();
        } else {
            std::cout << "Can not open binary sparse matrix file: " << filename << std::endl;
        }
    }

    ControlType
    getElectrodesControl(int n_wires, const NanowirePositionType &state1, const NanowirePositionType &state2) {
        ControlType control = 0;
        for (auto i = 0; i < n_wires; ++i) {
            if (abs(state1[i].first - state2[i].first) + abs(state1[i].second - state2[i].second) >= 2) {
                return 0;
            }

            if (state1[i] == state2[i]) {
                control |= (0b01 << 2 * (4 * state1[i].first + state1[i].second));
            } else {
                control |= (0b10 << 2 * (4 * state1[i].first + state1[i].second));
                control |= (0b01 << 2 * (4 * state2[i].first + state2[i].second));
            }
        }
        auto c = control;
        while (c != 0) {
            if ((c & 0b11) == 0b11)return 0;
            c = c >> 2;
        }
        return control;
    }

    void trimLeaf(NodePtr leaf){
        if(!leaf->getChildren().empty())return;
        auto parent = leaf->getParent();
        parent->removeChild(leaf);
        if(parent->getChildren().empty()){
            trimLeaf(parent);
        }
        leaf->setParent(nullptr);
        leaf.reset();
    }

    int getHeuristic(int n_wires,const NanowirePositionType& state1,const NanowirePositionType& state2){
        int value = 0;
        for(auto i=0;i<n_wires;++i) {
            auto temp_v = abs(state1[i].first - state2[i].first) + abs(state1[i].second - state2[i].second);
            value = std::max(value, temp_v);
        }
        return value;
    }

    int getHeuristic(int n_wires, IndexType index1, IndexType index2){

        int value = 0;
        for(auto i=0;i<n_wires;++i) {
            auto v1 = int(index1 & 0xf);
            auto v2 = int(index2 & 0xf);
            index1 = index1>>4;
            index2 = index2>>4;
            int temp_v = abs((v1>>2) - (v2>>2)) + abs((v1&0b11) - (v2&0b11));
            value = std::max(value, temp_v);
        }
        return value;
    }

    int getQuality(int n_wire,const NanowirePositionType &state)
    {
        int d1=1;
        int d2=0;
        for (auto i=0;i<n_wire-1;++i) {
            for(auto j=i+1;j<n_wire;++j) {
                auto t = (abs(state[i].first - state[j].first) + abs(state[i].second - state[j].second));
                d1 *= t;
                d2 += t;
            }
        }
        return -d1-d2;
    }

    int getQuality(int n_wire,IndexType state)
    {
        auto v = indexToElectrodeVector(n_wire,state);
        return getQuality(n_wire,v);
    }

    void destroyBranch(NodePtr node){
        if(node== nullptr)return;
        if(node->getParent()){
            node->getParent()->removeChild(node);
            node->setParent(nullptr);
        }
        auto children = node->getChildren();
        for(auto& child:children){
            destroyBranch(child);
        }
        node->getChildren().clear();
        node.reset();
    }


}
#endif //TRANSITIONMATRIXROUTE_UTILITY_HPP
