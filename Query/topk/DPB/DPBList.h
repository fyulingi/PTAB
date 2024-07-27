
#ifndef TOPK_DPB_ORDEREDLIST_H_
#define TOPK_DPB_ORDEREDLIST_H_
#include "../../../Util/Util.h"
#include "../CommonDP/Pool.h"
#include "../CommonDP/MinMaxHeap.hpp"
#include "../CommonDP/DynamicTrie.h"

class DPBList;
struct DpbElement{
  TYPE_ENTITY_LITERAL_ID node;
  unsigned int index;
  double cost;
  std::shared_ptr<DPBList> fq_;
  OnePointPredicatePtr predicates_;
  bool operator<(const DpbElement &other) const {return this->cost < other.cost;}
  bool operator!=(const DpbElement &other) const {
    return this->node!=other.node || this->cost!=other.cost || this->index !=other.index;
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
class DPBList{
 public:
  Pool pool_;
  DPBList()=default;
  // predicates_vec[i] is the predicates of child[i]
  // FQs and FRs will use this field.
  virtual OrderedListType Type(){return OrderedListType::UnDefined;};
  virtual void TryGetNext(unsigned int k)=0;
  virtual void GetResult(int i_th,std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr)=0;
  virtual ~DPBList(){};
};

class DPBFRIterator: public DPBList {
 private:
  // NodeOneChildVarPredicatesPtr type_predicates_;
  minmax::MinMaxHeap<DpbElement> queue_;
  std::vector<std::shared_ptr<DPBList>> pool_fqs_;
  std::vector<OnePointPredicatePtr> pool_predicates_;
  // std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DPBList>> fqs_map_;
 public:
  DPBFRIterator();
  virtual ~DPBFRIterator(){};
  OrderedListType Type() override {return OrderedListType::FR;};
  void TryGetNext(unsigned int k) override;

  // Insert a bulk of FQ iterator, all the same type
  void Insert(unsigned int k,
              TYPE_ENTITY_LITERAL_ID fq_id,
              std::shared_ptr<DPBList> FQ_iterator,
              OnePointPredicatePtr predicates_vec);

  static double DeltaCost(std::shared_ptr<DPBList>& node_pointer, int index);
  static bool NextEPoolElement(unsigned int k, std::shared_ptr<DPBList>& node_pointer, unsigned int index);
  virtual void GetResult(int i_th,std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
};

class DPBOWIterator: public DPBList{
 private:
  // sorting right after building, so don't need a heap to get top-k
  // the top-k result already in the pool
 public:
  DPBOWIterator() = default;
  virtual ~DPBOWIterator(){};
  OrderedListType Type() override {return OrderedListType::OW;};
  void TryGetNext(unsigned int k) override;

  // Insert a bulk of gStore Node ids and their scores
  void Insert(unsigned int k, const std::vector<TYPE_ENTITY_LITERAL_ID>& ids, const std::vector<double>& scores);
  void Insert(unsigned int k,const std::vector<TYPE_ENTITY_LITERAL_ID> &ids);
  virtual void GetResult(int i_th,std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
};

// FQ may have cost itself
class DPBFQIterator: public DPBList{
 private:
  double node_score_;
  minmax::MinMaxHeap<FqElement> queue_;
  std::vector<std::shared_ptr<DPBList>> fr_ow_iterators_;
  //std::vector<std::vector<unsigned int>> seq_list_;
  ePool e_pool_;
  DynamicTrie dynamic_trie_;
  std::vector<NodeOneChildVarPredicatesPtr> types_predicates_;
  TYPE_ENTITY_LITERAL_ID node_id_;
 public:

  virtual ~DPBFQIterator(){};
  inline void AddOneTypePredicate(NodeOneChildVarPredicatesPtr p){this->types_predicates_.push_back(std::move(p));}
  explicit DPBFQIterator(int k, TYPE_ENTITY_LITERAL_ID node_id, int child_type_num, double node_score):
      node_score_(node_score), dynamic_trie_(child_type_num,k),node_id_(node_id)
  {   this->fr_ow_iterators_.reserve(child_type_num);};
  void TryGetNext(unsigned int k) override;
  inline TYPE_ENTITY_LITERAL_ID GetNodeID(){return this->node_id_;};

  void Insert(std::vector<std::shared_ptr<DPBList>> FR_OW_iterators);
  void Insert(std::shared_ptr<DPBList> FR_OW_iterator);
  static bool NextEPoolElement(unsigned int k, std::shared_ptr<DPBList>& node_pointer, unsigned int index);
  static double DeltaCost(std::shared_ptr<DPBList>& FR_OW_iterator, int index);
  virtual void GetResult(int i_th,std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
  virtual OrderedListType Type() override {return OrderedListType::FQ;};
};

#endif //TOPK_DPB_ORDEREDLIST_H_
