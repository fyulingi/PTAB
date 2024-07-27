
// 需要进行的两个优化分别是：
// 1. explore 的时候 explore 到一个即停止 done
//    可能要重新思考 uncertainty level的写法了 done
// 2. try get next 的时候引入 second lower bound 的概念 [放弃这个]

#ifndef PRIVATETOPK_QUERY_TOPK_Dfs_PLIST_H_
#define PRIVATETOPK_QUERY_TOPK_Dfs_PLIST_H_

#include "../../../Util/Util.h"
#include "../CommonDP/Pool.h"
#include "../CommonDP/MinMaxHeap.hpp"
#include "./CompressedVector.h"
#include "./HopIndex.h"
#include "../../../Query/IDList.h"
#include "../TopKSearchPlan.h"
#include "../TopKUtil.h"

class DfsList;

struct NodeCaches{
  std::unique_ptr<std::unordered_map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DfsList>>[]> var_valid_ids;
  std::unique_ptr<std::unordered_set<TYPE_ENTITY_LITERAL_ID>[]> var_invalid_ids;
};

struct DfsHelpInfo{
  NodeCaches *node_caches;
  HopIndex* hop_index;
  TopKTreeNode** node_plan;
  TopKUtil::Env *env;
  // below is the contents of TopKUtil::Env
  KVstore *kv_store;
  BasicQuery *basic_query;
  shared_ptr<BGPQuery> bgp_query;
  shared_ptr<unordered_map<TYPE_ENTITY_LITERAL_ID,shared_ptr<IDList>>> id_caches;
  int k;
  // store the non-tree edges to speed up result enumeration
  std::shared_ptr<std::vector<std::map<TYPE_ENTITY_LITERAL_ID,std::set<TYPE_ENTITY_LITERAL_ID> >>> non_tree_edges_lists_;
  std::shared_ptr<unordered_map<std::string, double>> coefficients;
  shared_ptr<Transaction> txn;
  std::shared_ptr<stringstream> ss;
  std::vector<std::vector<bool>> ow_coefficient_counts_;
  std::vector<bool> subtree_has_coefficient_;
  std::unique_ptr<size_t[]> max_leaf_distance_;
  bool DfsNew;
};

enum class DfsValue{Exact,Range};
struct DfsElement{
  TYPE_ENTITY_LITERAL_ID node;
  unsigned int index;
  DfsValue value_type;
  double cost;
  bool operator<(const DfsElement &other) const {
    if(this->cost < other.cost)
      return true;
    if(this->cost == other.cost)
      return this->value_type==DfsValue::Exact;
    return false;
  }
  bool operator!=(const DfsElement &other) const {
    return this->node!=other.node || this->value_type!= other.value_type|| this->cost!=other.cost || this->index !=other.index;
  }
};

struct DfsFqElement{
 public:
  sequence seq;
  DfsValue value_type;
  double cost;
  DfsFqElement() = default;
  DfsFqElement(sequence s,DfsValue v, double c):seq(std::move(s)),value_type(v),cost(c){}
  bool operator<(const DfsFqElement &other) const {
    if(this->cost < other.cost)
      return true;
    if(this->cost == other.cost)
      return this->value_type==DfsValue::Exact;
    return false;
  }

  bool operator!=(const DfsFqElement& other) const {
    return this->cost != other.cost ||this->value_type!= other.value_type ||!std::equal(this->seq.begin(),this->seq.end(),other.seq.begin());
  }
};
/**
 * Structure
 * FR
 *   FQ
 *   predicates with this FQ
 *
 * FQ [OW first and FRs last]
 *   OWs
 *   predicates information with all nodes in OW
 *   FRs
 *   no predicates information
 * OW
 *   no predicates information
 */

/* This structure is designed for
  selecting top-k min element */
class DfsList{
 public:
  bool exhausted_;
  bool not_fully_explored_;
  bool has_estimated_;
  double estimate_[D_HOP_LOG];
  Pool pool_;

  virtual bool Empty() = 0;

  // it should let the 'pool_' counter increase 1 , making one more exact match
  // and it will update the estimate, as the estimate only reflect the fact
  // about pool[n], but n has increased just
  virtual void TryGetNext(DfsHelpInfo *dfs_help_info)=0;

  DfsList(): has_estimated_(false), exhausted_(false),not_fully_explored_(true){};
  // predicates_vec[i] is the predicates of child[i]
  // FQs and FRs will use this field.
  virtual OrderedListType Type(){return OrderedListType::UnDefined;};
  virtual void GetResult(int i_th,std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr)=0;
  virtual ~DfsList(){};
  /**
   * Estimate() the lower bound of top-1 answer.
   */
  virtual void Estimate(){};
  virtual int GetVarId(){ return -1;};
  /**
   * @return whether the Iterator will return no more values
   */
  inline bool Exhausted() {return this->exhausted_; }
  inline void SetExhausted(bool v){ this->exhausted_ = v; }
};


#endif //PRIVATETOPK_QUERY_TOPK_Dfs_PLIST_H_
