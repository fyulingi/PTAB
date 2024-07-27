
#ifndef PRIVATETOPK_QUERY_TOPK_Dfs_DfsFQITERATOR_H_
#define PRIVATETOPK_QUERY_TOPK_Dfs_DfsFQITERATOR_H_
#include "./DfsList.h"
#include "./HopIndex.h"

class DfsFRIterator;
// FQ may have cost itself
class DfsFQIterator: public DfsList{
 private:
  double node_score_;
  minmax::MinMaxHeap<DfsFqElement> queue_;
  // range score
  std::vector<std::shared_ptr<DfsList>> fr_ow_iterators_;
  // exact score and it is strictly refers to one pool_
  ePool e_pool_;
  CompressedVector compressed_vector_;
  // DynamicTrie dynamic_trie_;
  std::vector<NodeOneChildVarPredicatesPtr> types_predicates_;
  TYPE_ENTITY_LITERAL_ID node_id_;
  int var_id_;
 public:
  bool Empty() {
    if(!this->pool_.empty())
      return false;
    return  this->fr_ow_iterators_.empty();
  };
  virtual ~DfsFQIterator(){};
  inline std::shared_ptr<DfsList> GetChild(unsigned int i_th){return this->fr_ow_iterators_[i_th];}
  inline void AddOneTypePredicate(NodeOneChildVarPredicatesPtr p){this->types_predicates_.push_back(std::move(p));}
  explicit DfsFQIterator(int k, TYPE_ENTITY_LITERAL_ID node_id, int var_id, int child_type_num, double node_score):
      node_score_(node_score),
      // dynamic_trie_(child_type_num,k),
      compressed_vector_(k,child_type_num),
      node_id_(node_id),var_id_(var_id)
  {
    this->fr_ow_iterators_.reserve(child_type_num);
  };
  void TryGetNext(DfsHelpInfo *dfs_help_info);
  inline TYPE_ENTITY_LITERAL_ID GetNodeID(){return this->node_id_;};
  int GetVarId() override {return this->var_id_;}
  void Insert(std::vector<std::shared_ptr<DfsList>> FR_OW_iterators);
  void Insert(std::shared_ptr<DfsList> FR_OW_iterator);
  virtual void GetResult(int i_th,std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
  void Estimate() override;
  virtual OrderedListType Type() override {return OrderedListType::FQ;};
};


#endif //PRIVATETOPK_QUERY_TOPK_Dfs_DfsFQITERATOR_H_
