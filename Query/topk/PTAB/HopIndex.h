
#ifndef PRIVATETOPK_QUERY_TOPK_DFS_HOPINDEX_H_
#define PRIVATETOPK_QUERY_TOPK_DFS_HOPINDEX_H_
#include "../../../Util/Util.h"
#include "../../../KVstore/KVstore.h"
#include "../CommonDP/MinMaxHeap.hpp"
#include <iostream>
#define HOP_INDEX_DEBUG_INFO

//500 GB
// 200 GB for BSBM
constexpr size_t memory_avaiable = 250LL * 1000LL * 1000LL * 1000LL;
/**
 * Predicate file should be like:
 * [line 0-n-1] predicate 0-n-1
 */
class TopTwoInfo{
 public:
  // the value of the minimum ,second-minimum, maximum, second-maximum
  double min_1,min_2,max_1,max_2;
  TYPE_ENTITY_LITERAL_ID min_1_node, min_2_node, max_1_node , max_2_node;
  TopTwoInfo(){
    min_1_node=INVALID_ENTITY_LITERAL_ID;
    min_2_node = INVALID_ENTITY_LITERAL_ID;
    max_1_node= INVALID_ENTITY_LITERAL_ID;
    max_2_node = INVALID_ENTITY_LITERAL_ID;
  };
  TopTwoInfo(TYPE_ENTITY_LITERAL_ID initial_id ,double min_v, double max_v)
  {
    min_1_node= initial_id;
    min_1 = min_v;
    min_2_node = INVALID_ENTITY_LITERAL_ID;
    max_1_node= initial_id;
    max_1 = max_v;
    max_2_node = INVALID_ENTITY_LITERAL_ID;
  }
  void AddMinNode(TYPE_ENTITY_LITERAL_ID id, double min_v);
  void AddMaxNode(TYPE_ENTITY_LITERAL_ID id, double max_v);
  void AddAnotherNode(const TopTwoInfo &other);
  std::tuple<bool,double,bool,double> GetRange(TYPE_ENTITY_LITERAL_ID id) const;
};

enum class RangeState{Exact,UnExact,NoExist};
char RangeStateToChar(RangeState state);
RangeState CharToRangeState(char x);

struct NodeHopInfo{
  unsigned int predicate_num_;
  std::unique_ptr<TYPE_PREDICATE_ID[]> predicate_ids;
  // corresponding_nodes[i] records the literal id of the object of predicate[i]
  // if only one literal exists
  std::unique_ptr<TYPE_ENTITY_LITERAL_ID []> corresponding_nodes;

  // for predicate_ids[i]
  // min_values_[D_HOP*i] -> min_values_[D_HOP*i + log(2)(D_HOP) - 1]
  // refers to the min value of the predicate in 1,2,4,...,D_HOP
  std::unique_ptr<double[]> min_values_;
  std::unique_ptr<double[]> max_values_;
  // mask_[D_HOP_LOG*i + j] means
  // j = 0 ,mask_[D_HOP_LOG*i + j] == true means the predicate has only one literal
  // j >= 1 min_values_[D_HOP_LOG * i + j ] and max_values_[D_HOP_LOG * i + j ] is used
  std::unique_ptr<std::vector<RangeState>> mask_;
  inline unsigned int GetPrePosition(TYPE_PREDICATE_ID pre_id) const
  {
    auto lower = std::lower_bound(predicate_ids.get(), predicate_ids.get() + this->predicate_num_, pre_id);
    auto pos = lower-predicate_ids.get();
    if(pos == this->predicate_num_ || *lower != pre_id)
      return -1;
    return pos;
  }
};

/**
 * The File is organized as follows:
 * [0- meta_data_offset ]  entity_num_,predicate_num_,type_pre_id_
 * data of entities_type_
 * data of entity_has_type
 * data of predicates_to_number_
 * data of hop_info_
 */

class HopIndex {
  bool loaded;
  size_t valid_num;
  TYPE_ENTITY_LITERAL_ID entity_num_;
  TYPE_PREDICATE_ID predicate_num_;
  // the ID corresponding to predicate like, <http://www.w3.org/1999/02/22-rdf-syntax-ns#type>
  TYPE_PREDICATE_ID type_pre_id_;
  static const  TYPE_ENTITY_LITERAL_ID TypeInvalid = -1;

  std::unique_ptr<TYPE_ENTITY_LITERAL_ID[]> entities_type_;
  // predicates whose object has a type that can be converted into numbers
  std::set<TYPE_PREDICATE_ID> predicates_to_number_;
  // std::unique_ptr<NodeHopInfo[]> hop_info_;
  size_t predicate_total;

  // return the max/min value of the set of all predicates in ids
  std::unique_ptr<std::map<TYPE_PREDICATE_ID,TopTwoInfo>> MergeSameTypeNodes(const std::set<TYPE_ENTITY_LITERAL_ID>& ids)const;
  size_t CountAttrNum(const std::set<TYPE_ENTITY_LITERAL_ID>& ids)const;
 public:
  // One node will occupy
  //     content            numbers
  //   predicate_ids        predicate_num_
  // corresponding_nodes    predicate_num_
  //    min_values_         predicate_num_ * D_HOP_LOG
  //    max_values_         predicate_num_ * D_HOP_LOG
  //    mask_               predicate_num_ * D_HOP_LOG
  std::unique_ptr<size_t[]> entity_offset; // entity_num_
  std::unique_ptr<size_t[]> predicates_num;// entity_num_
  std::unique_ptr<TYPE_PREDICATE_ID[]> predicate_ids;// \sum predicates_num
  std::unique_ptr<TYPE_ENTITY_LITERAL_ID []> corresponding_nodes;// \sum predicates_num
  std::unique_ptr<double[]> min_values_;// \sum predicates_num * D_HOP_LOG
  std::unique_ptr<double[]> max_values_;// \sum predicates_num * D_HOP_LOG
  std::unique_ptr<RangeState[]> mask_;// \sum predicates_num * D_HOP_LOG

  HopIndex(TYPE_ENTITY_LITERAL_ID entity_num,TYPE_PREDICATE_ID predicate_num):
      loaded(false),entity_num_(entity_num),predicate_num_(predicate_num){
    cout<<"Entity num:"<<entity_num<<" predicates num:"<<predicate_num<<endl;
    this->entities_type_ = std::unique_ptr<TYPE_ENTITY_LITERAL_ID[]> (new TYPE_ENTITY_LITERAL_ID[entity_num]);
    memset( this->entities_type_.get(),-1,entity_num*sizeof(TYPE_ENTITY_LITERAL_ID) );
  };
  void LoadPredicates(KVstore *kv_store,std::string file_name);
  void BuildNodeInfo(KVstore *kv_store);
  void BuildHopInfo(KVstore *kv_store);
  void Build(KVstore *kv_store);
//  std::unique_ptr<NodeHopInfo[]> & GetHopIndexRef(){return this->hop_info_;}
//  inline const NodeHopInfo& GetNodeInfo(TYPE_ENTITY_LITERAL_ID id){
//    return this->hop_info_[id];
//  }
  inline bool IsValidPredicateToNumber(TYPE_PREDICATE_ID pre) {
    return this->predicates_to_number_.find(pre) != this->predicates_to_number_.end();
  }
  void Save(std::string file_name);
  void Load(std::string file_name);
  size_t GetOffset(TYPE_ENTITY_LITERAL_ID node_id){
    return this->entity_offset[node_id];
  }
  size_t GetOffset(TYPE_ENTITY_LITERAL_ID node_id,TYPE_PREDICATE_ID pre_id)
  {
    auto node_offset = GetOffset(node_id);
    auto pre_num = this->predicates_num[node_id];
    if(pre_num == 0) return -1;
    auto p_start = this->predicate_ids.get() + node_offset;
    auto p_end = p_start + pre_num;
    auto lower = std::lower_bound(p_start, p_end, pre_id);
    if(lower == p_end)
      return -1;
    auto pos = lower - this->predicate_ids.get();
    if(pos == this->predicate_total ||*lower != pre_id)
      return -1;
    return pos;
  }
};

#endif //PRIVATETOPK_QUERY_TOPK_DFS_HOPINDEX_H_
