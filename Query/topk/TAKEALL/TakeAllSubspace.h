
#ifndef TAKE_ALL_SUBSPACE_H_
#define TAKE_ALL_SUBSPACE_H_
#include "../../../Util/Util.h"
#include "../CommonDP/Pool.h"
#include "../CommonDP/MinMaxHeap.hpp"

class TakeAllSubspace;
class TakeAllSubspaceLess;
class TakeAllFQIterator;
using TakeAllSpaceHeap = minmax::MinMaxHeap<std::unique_ptr<TakeAllSubspace>,std::vector<std::unique_ptr<TakeAllSubspace>>,TakeAllSubspaceLess>;

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
  selecting top-k min TakeAllElement */
class TakeAllList: public std::enable_shared_from_this<TakeAllList>{
 public:
  std::unique_ptr<TakeAllElement> top_element_ptr_;
  TakeAllList()=default;
  // predicates_vec[i] is the predicates of child[i]
  // FQs and FRs will use this field.
  virtual OrderedListType Type(){return OrderedListType::UnDefined;};
  virtual void GetFirst(unsigned int k)=0;
  virtual void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr)=0;
  virtual ~TakeAllList()=default;

  // cut a node from this iterator and this TakeAllList will changed to two new List
  // the origin will remain
  virtual void
  SplitAllIntoSpaces(unsigned int k,
                     double space_best_score,
                     std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                     std::deque<std::shared_ptr<TakeAllFQIterator>> &fq_stack,
                     std::deque<unsigned int> &iterate_times,
                     TakeAllSpaceHeap *collection)=0;
  virtual void
  MergeUpSubSpace(unsigned int k,
                    double space_best_score,
                    std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                    std::deque<std::shared_ptr<TakeAllFQIterator>> &fq_stack,
                    std::deque<unsigned int> &iterate_times,
                    TakeAllSpaceHeap *collection)=0;

 protected:
  template <typename Derived>
  std::shared_ptr<Derived> shared_from_base()
  {
    return std::static_pointer_cast<Derived>(shared_from_this());
  }
};

struct TakeAllfrElement{
  TYPE_ENTITY_LITERAL_ID node;
  double cost;
  std::shared_ptr<TakeAllList> fq_;
  OnePointPredicatePtr predicate_;
  bool operator<(const TakeAllfrElement &other) const {return this->cost < other.cost;}
  bool operator!=(const TakeAllfrElement &other) const {
    return this->node!=other.node || this->cost!=other.cost;
  }
};

class TakeAllFRIterator: public TakeAllList {
 private:

  std::vector<TYPE_ENTITY_LITERAL_ID> fq_ids_;
  std::vector<OnePointPredicatePtr> predicates_;
  // NodeOneChildVarPredicatesPtr type_predicates_;
  std::vector<std::shared_ptr<TakeAllList>> fqs_;
  // std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<TakeAllList>> fqs_map_;
 public:
  std::unique_ptr<TakeAllfrElement> top_element_ptr_;
  TakeAllFRIterator();
  ~TakeAllFRIterator(){};
  OrderedListType Type() override {return OrderedListType::FR;};
  void GetFirst(unsigned int k) override;

  // Insert a bulk of FQ iterator, all the same type
  void Insert(unsigned int k,
              TYPE_ENTITY_LITERAL_ID fq_id,
              std::shared_ptr<TakeAllList> FQ_iterator,
              OnePointPredicatePtr predicates_vec);

  void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                 NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;

  void
  SplitAllIntoSpaces(unsigned int k,
                     double space_best_score,
                     std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                     std::deque<std::shared_ptr<TakeAllFQIterator>> &fq_stack,
                     std::deque<unsigned int> &iterate_times,
                     TakeAllSpaceHeap *collection) override;

  // split fr which has only one node
  void SplitTopFq(unsigned int k,
                  double space_best_score,
                  std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                  std::deque<std::shared_ptr<TakeAllFQIterator>> &fq_stack,
                  std::deque<unsigned int> &iterate_times,
                  TakeAllSpaceHeap *collection);

  void
  MergeUpSubSpace(unsigned int k,
                  double space_best_score,
                  std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                  std::deque<std::shared_ptr<TakeAllFQIterator>> &fq_stack,
                  std::deque<unsigned int> &iterate_times,
                  TakeAllSpaceHeap *collection) override;
};

class TakeAllOWIterator: public TakeAllList{
 private:
  std::vector<TakeAllElement> pool_;
 public:
  TakeAllOWIterator() = default;
  ~TakeAllOWIterator(){};
  OrderedListType Type() override {return OrderedListType::OW;};
  void GetFirst(unsigned int k) override;

  // Insert a bulk of gStore Node ids and their scores
  void Insert(unsigned int k, const std::vector<TYPE_ENTITY_LITERAL_ID>& ids, const std::vector<double>& scores);
  void Insert(unsigned int k,const std::vector<TYPE_ENTITY_LITERAL_ID> &ids);
  void RankedInsert(unsigned int k, const std::vector<TYPE_ENTITY_LITERAL_ID>& ids, const std::vector<double>& scores);
  void RankedInsert(unsigned int k,TYPE_ENTITY_LITERAL_ID id,double score);
  void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
  void
  SplitAllIntoSpaces(unsigned int k,
                     double space_best_score,
                     std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                     std::deque<std::shared_ptr<TakeAllFQIterator>> &fq_stack,
                     std::deque<unsigned int> &iterate_times,
                     TakeAllSpaceHeap *collection) override;

  void
  MergeUpSubSpace(unsigned int k,
                  double space_best_score,
                  std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                  std::deque<std::shared_ptr<TakeAllFQIterator>> &fq_stack,
                  std::deque<unsigned int> &iterate_times,
                  TakeAllSpaceHeap *collection)override;
};

// FQ may have cost itself
class TakeAllFQIterator: public TakeAllList{
 private:
  double node_score_;
  std::vector<std::shared_ptr<TakeAllList>> fr_ow_iterators_;
  std::vector<NodeOneChildVarPredicatesPtr> types_predicates_;
  TYPE_ENTITY_LITERAL_ID node_id_;
 public:

  virtual ~TakeAllFQIterator(){};
  void AddOneTypePredicate(NodeOneChildVarPredicatesPtr p){this->types_predicates_.push_back(p);}
  explicit TakeAllFQIterator(int k, TYPE_ENTITY_LITERAL_ID node_id, int child_type_num, double node_score):
      node_score_(node_score),node_id_(node_id)
  {   this->fr_ow_iterators_.reserve(child_type_num);};
  void GetFirst(unsigned int k) override;
  inline TYPE_ENTITY_LITERAL_ID GetNodeID(){return this->node_id_;};

  void Insert(std::vector<std::shared_ptr<TakeAllList>> FR_OW_iterators);
  void Insert(std::shared_ptr<TakeAllList> FR_OW_iterator);
  inline unsigned int ChildTypeNum() {return this->fr_ow_iterators_.size();}
  virtual void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
  virtual OrderedListType Type() override {return OrderedListType::FQ;};
  void
  SplitAllIntoSpaces(unsigned int k,
                     double space_best_score,
                     std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                     std::deque<std::shared_ptr<TakeAllFQIterator>> &fq_stack,
                     std::deque<unsigned int> &iterate_times,
                     TakeAllSpaceHeap *collection) override;
  void
  MergeUpSubSpace(unsigned int k,
                  double space_best_score,
                  std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                  std::deque<std::shared_ptr<TakeAllFQIterator>> &fq_stack,
                  std::deque<unsigned int> &iterate_times,
                  TakeAllSpaceHeap *collection)override;
};

class TakeAllSubspace{
 public:
  std::vector<TYPE_ENTITY_LITERAL_ID> one_node_space_;
  std::deque<std::shared_ptr<TakeAllFQIterator>> fq_stack_;
  // records which i-th child of fq_stack[i] is being visited
  std::deque<unsigned int> iterate_times_;
  // end_node_ must be FR/OW
  std::shared_ptr<TakeAllList> end_node_;
  double best_score_;
  bool operator<(const TakeAllSubspace& other) const
  {
    return this->best_score_ < other.best_score_;
  }
};

class TakeAllSubspaceLess{
 public:
  inline bool operator()(const std::unique_ptr<TakeAllSubspace>& p1,const std::unique_ptr<TakeAllSubspace>& p2) const
  {
    return p1->best_score_<p2->best_score_;
  }
};


#endif //TAKE_ALL_SUBSPACE_H_
