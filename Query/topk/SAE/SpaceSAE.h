
#ifndef TOPK_SpaceSAE_H_
#define TOPK_SpaceSAE_H_
#include "../../../Util/Util.h"
#include "../CommonDP/Pool.h"
#include "../CommonDP/MinMaxHeap.hpp"

class SpaceSAE;
class SpaceSAELess;
class SAEFQIterator;
using SpaceHeapSAE = minmax::MinMaxHeap<std::unique_ptr<SpaceSAE>,std::vector<std::unique_ptr<SpaceSAE>>,SpaceSAELess>;

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
  selecting top-k min SAEElement */
class SAEList{
 public:
  std::unique_ptr<SAEElement> top_element_ptr_;
  SAEList()=default;
  // predicates_vec[i] is the predicates of child[i]
  // FQs and FRs will use this field.
  virtual OrderedListType Type(){return OrderedListType::UnDefined;};
  virtual void GetFirst(unsigned int k)=0;
  virtual void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr)=0;
  virtual ~SAEList()=default;

  // 从上往下深度优先遍历的过程中，切割子空间的过程
  // cut a node from this iterator and this SAEList will changed to two new List
  // the origin will remain
  virtual void
  SplitIntoTwoSpace(unsigned int k,
                    double space_best_score,
                    std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                    std::unique_ptr<std::stack<std::shared_ptr<SAEList>>> fr_ow_iterators,
                    unsigned int iterate_time,
                    SpaceHeapSAE *collection)=0;
};

struct SAEfrElement{
  TYPE_ENTITY_LITERAL_ID node;
  double cost;
  std::shared_ptr<SAEList> fq_;
  OnePointPredicatePtr predicate_;
  bool operator<(const SAEfrElement &other) const {return this->cost < other.cost;}
  bool operator!=(const SAEfrElement &other) const {
    return this->node!=other.node || this->cost!=other.cost;
  }
};

class SAEFRIterator: public SAEList {
 private:
  std::vector<TYPE_ENTITY_LITERAL_ID> fq_ids_;
  std::vector<OnePointPredicatePtr> predicates_;
  std::vector<std::shared_ptr<SAEList>> fqs_;
  std::shared_ptr<SAEFRIterator> next_;
 public:
  std::unique_ptr<SAEfrElement> top_element_ptr_;
  SAEFRIterator();
  ~SAEFRIterator(){};
  OrderedListType Type() override {return OrderedListType::FR;};
  void GetFirst(unsigned int k) override;

  // Insert a bulk of FQ iterator, all the same type
  void Insert(unsigned int k,
              TYPE_ENTITY_LITERAL_ID fq_id,
              std::shared_ptr<SAEList> FQ_iterator,
              OnePointPredicatePtr predicates_vec);

  void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                 NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;

  void
  SplitIntoTwoSpace(unsigned int k,
                    double space_best_score,
                    std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                    std::unique_ptr<std::stack<std::shared_ptr<SAEList>>> fr_ow_iterators,
                    unsigned int iterate_time,
                    SpaceHeapSAE *collection) override;

  // split fr which has only one node
  void SplitTopFq(unsigned int k,
                  double space_best_score,
                  std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                  SpaceHeapSAE *collection);

};

class SAEOWIterator: public SAEList{
 private:
  std::vector<SAEElement> pool_;
  std::shared_ptr<SAEOWIterator> next_;
 public:
  SAEOWIterator() = default;
  ~SAEOWIterator(){};
  OrderedListType Type() override {return OrderedListType::OW;};
  void GetFirst(unsigned int k) override;

  // Insert a bulk of gStore Node ids and their scores
  void Insert(unsigned int k, const std::vector<TYPE_ENTITY_LITERAL_ID>& ids, const std::vector<double>& scores);
  void Insert(unsigned int k,const std::vector<TYPE_ENTITY_LITERAL_ID> &ids);
  void RankedInsert(unsigned int k, const std::vector<TYPE_ENTITY_LITERAL_ID>& ids, const std::vector<double>& scores);
  void RankedInsert(unsigned int k, TYPE_ENTITY_LITERAL_ID id,double score);
  void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
  void
  SplitIntoTwoSpace(unsigned int k,
                    double space_best_score,
                    std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                    std::unique_ptr<std::stack<std::shared_ptr<SAEList>>> fr_ow_iterators,
                    unsigned int iterate_time,
                    SpaceHeapSAE *collection) override;

};

// SpaceSAE must have 2 ability：
// 1. acquire the best solution
// 2. split to more sub-spaces
// In SAE algorithm, Space only concerns about the fr/ow this layer
class SpaceSAE{
 public:
  std::vector<TYPE_ENTITY_LITERAL_ID> one_node_space_;
  // iterated_items_ means fq_ has iterated 'iterated_items_' times
  // if fq has 3 children R1 R2 R3, iterated_items_ = 1
  // it means the space is {r1} x R2 x R3
  // r1(may have multiple nodes, say |r1| = n ) stores in one_node_space_[0...(n-1)]
  unsigned int iterated_items_;
  // R2 R3 store in stack bottom [R3 R2] top
  std::unique_ptr<std::stack<std::shared_ptr<SAEList>>> fr_ow_iterators_;
  double best_score_;
  SpaceSAE(std::vector<TYPE_ENTITY_LITERAL_ID> one_node_space,
           unsigned int iterated_items,
           std::unique_ptr<std::stack<std::shared_ptr<SAEList>>> fr_ow_iterators,
           double best_score):one_node_space_(std::move(one_node_space)),
                              iterated_items_(iterated_items),
                              fr_ow_iterators_(std::move(fr_ow_iterators)),
                              best_score_(best_score)
           {}
  void Set(std::vector<TYPE_ENTITY_LITERAL_ID> one_node_space,
           unsigned int iterated_items,
           std::unique_ptr<std::stack<std::shared_ptr<SAEList>>>  fr_ow_iterators,
           double best_score);
  SpaceSAE(SpaceSAE &&other):one_node_space_ (std::move(other.one_node_space_)),
                             iterated_items_(std::move(other.iterated_items_)),
                             fr_ow_iterators_(std::move(other.fr_ow_iterators_)),
                             best_score_(std::move(other.best_score_))
  {};
  bool operator<(const SpaceSAE& other) const
  {
    return this->best_score_ < other.best_score_;
  }
};

class SpaceSAELess{
 public:
  inline bool operator()(const std::unique_ptr<SpaceSAE>& p1,const std::unique_ptr<SpaceSAE>& p2) const
  {
    return p1->best_score_<p2->best_score_;
  }
};

// FQ may have cost itself
 class SAEFQIterator: public SAEList{
 private:
  double node_score_;
  std::vector<std::shared_ptr<SAEList>> fr_ow_iterators_;
  std::vector<NodeOneChildVarPredicatesPtr> types_predicates_;
  TYPE_ENTITY_LITERAL_ID node_id_;
  SpaceHeapSAE* small_queue_;
 public:
  virtual ~SAEFQIterator(){
    if(small_queue_!= nullptr)
      delete small_queue_;
      small_queue_ = nullptr;
  };
  void AddOneTypePredicate(NodeOneChildVarPredicatesPtr p){this->types_predicates_.push_back(p);}
  explicit SAEFQIterator(int k, TYPE_ENTITY_LITERAL_ID node_id, int child_type_num, double node_score):
      node_score_(node_score),node_id_(node_id),small_queue_(nullptr)
  {   this->fr_ow_iterators_.reserve(child_type_num);};
  void GetFirst(unsigned int k) override;
  inline TYPE_ENTITY_LITERAL_ID GetNodeID(){return this->node_id_;};

  void Insert(std::vector<std::shared_ptr<SAEList>> FR_OW_iterators);
  void Insert(std::shared_ptr<SAEList> FR_OW_iterator);
  inline unsigned int ChildTypeNum() {return this->fr_ow_iterators_.size();}
  virtual void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
  virtual OrderedListType Type() override {return OrderedListType::FQ;};
  // 谓词这里的ID还没搞好
  void
  SplitIntoTwoSpace(unsigned int k,
                    double space_best_score,
                    std::vector<TYPE_ENTITY_LITERAL_ID>& record_up_layer,
                    std::unique_ptr<std::stack<std::shared_ptr<SAEList>>> fr_ow_iterators,
                    unsigned int iterate_time,
                    SpaceHeapSAE *collection) override;

};


#endif //TOPK_SpaceSAE_H_
