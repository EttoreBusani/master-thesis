#ifndef PERMUTATIONS_HPP
#define PERMUTATIONS_HPP

#include<vector>

// This file contains the functions used to generate a vector
// of all the permutations of the names of the platform applications.

// Append the input-value to each element of the input-vector.
// INPUTS: el -> element to be appended
//         vec -> vector of vectors to be modified.
// OUTPUTS: -
template<typename T>
void appendEl(T el, std::vector<std::vector<T>>& vec){
  for(auto it=vec.begin(); it!=vec.cend(); ++it){
    it->push_back(el);
  }
}

// Appends the 2° input-vector to the 1° input-vector.
// INPUTS: vec1 -> vector to be modified
//         vec2 -> vector to be appended.
template<typename T>
void append(std::vector<T>& vec1, const std::vector<T>& vec2){
  for(auto it=vec2.cbegin(); it!=vec2.cend(); ++it) vec1.push_back(*it);
}

// Returns the vector containing all the permutations of the input-vector
// (Recoursive function)
// INPUTS: vec -> vector of elements to be permuted.
// OUTPUTS vector containing permutations.
template<typename T>
std::vector<std::vector<T>> getPermutations(const std::vector<T>& vec){
  // Base cases.
  if(vec.size() <= 1) return {vec};
  else if(vec.size() == 2) return {vec, {vec[1], vec[0]}};
  // Recoursion
  else{
    std::vector<std::vector<T>> res = {};
    for(auto el : vec) {
      std::vector<T> recVec = {};
      for(auto el1 : vec){
        if(el1!=el) recVec.push_back(el1);
      }
      auto curr = getPermutations(recVec);
      appendEl(el,curr);
      append(res,curr);
    }
    return res;
  }
}

// Returns the vector containing the permutations of the first n
// elements of input vector:
// example: getPermutations({1,2,3},1) returns {{1,2,3},{2,1,3},{3,1,2}}
// (Recoursive function)
// INPUTS: vec -> vector of elements to be permuted.
//          n -> n. of elements to be permuted
// OUTPUTS vector containing permutations.
template<typename T>
std::vector<std::vector<T>> getPermutations(const std::vector<T> vec, int n){
  if(vec.size()==n) return getPermutations(vec);
  else{
    std::vector<std::vector<T>> res={};
    if(n==0) return{vec};
    else{
      for(auto el : vec){
        std::vector<T> newVec={};
        for(auto el1:vec) if(el1!=el) newVec.push_back(el1);
        auto res1 = getPermutations(newVec,n-1);
        for(auto el2 : res1){
          std::vector<T> perm={el};
          for(auto el3 : el2) perm.push_back(el3);
          res.push_back(perm);
        }

      }
    }
    return res;
  }
}


#endif
