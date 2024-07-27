
#ifndef TOPK_SUBSPACE_H_
#define TOPK_SUBSPACE_H_
#include "../../../Util/Util.h"
#include "../CommonDP/Pool.h"
#include "../CommonDP/MinMaxHeap.hpp"

class Subspace;
class SubspaceLess;
class kTPMFQIterator;
using SpaceHeap = minmax::MinMaxHeap<std::unique_ptr<Subspace>,std::vector<std::unique_ptr<Subspace>>,SubspaceLess>;

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
  selecting top-k min kTPMElement */
class kTPMList: public std::enable_shared_from_this<kTPMList>{
 public:
  std::unique_ptr<kTPMElement> top_element_ptr_;
  kTPMList()=default;
  // predicates_vec[i] is the predicates of child[i]
  // FQs and FRs will use this field.
  virtual OrderedListType Type(){return OrderedListType::UnDefined;};
  virtual void GetFirst(unsigned int k)=0;
  virtual void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr)=0;
  virtual ~kTPMList()=default;

  // cut a node from this iterator and this kTPMList will changed to two new List
  // the origin will remain
  virtual void
  SplitIntoTwoSpace(unsigned int k,
                    double space_best_score,
                    std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                    std::deque<std::shared_ptr<kTPMFQIterator>> &fq_stack,
                    std::deque<unsigned int> &iterate_times,
                    SpaceHeap *collection)=0;
  virtual void
  MergeUpSubSpace(unsigned int k,
                    double space_best_score,
                    std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                    std::deque<std::shared_ptr<kTPMFQIterator>> &fq_stack,
                    std::deque<unsigned int> &iterate_times,
                    SpaceHeap *collection)=0;

 protected:
  template <typename Derived>
  std::shared_ptr<Derived> shared_from_base()
  {
    return std::static_pointer_cast<Derived>(shared_from_this());
  }
};

struct kTPMfrElement{
  TYPE_ENTITY_LITERAL_ID node;
  double cost;
  std::shared_ptr<kTPMList> fq_;
  OnePointPredicatePtr predicate_;
  bool operator<(const kTPMfrElement &other) const {return this->cost < other.cost;}
  bool operator!=(const kTPMfrElement &other) const {
    return this->node!=other.node || this->cost!=other.cost;
  }
};

class kTPMFRIterator: public kTPMList {
 private:

  std::vector<TYPE_ENTITY_LITERAL_ID> fq_ids_;
  std::vector<OnePointPredicatePtr> predicates_;
  // NodeOneChildVarPredicatesPtr type_predicates_;
  std::vector<std::shared_ptr<kTPMList>> fqs_;
  // std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<kTPMList>> fqs_map_;
 public:
  std::unique_ptr<kTPMfrElement> top_element_ptr_;
  kTPMFRIterator();
  ~kTPMFRIterator(){};
  OrderedListType Type() override {return OrderedListType::FR;};
  void GetFirst(unsigned int k) override;

  // Insert a bulk of FQ iterator, all the same type
  void Insert(unsigned int k,
              TYPE_ENTITY_LITERAL_ID fq_id,
              std::shared_ptr<kTPMList> FQ_iterator,
              OnePointPredicatePtr predicates_vec);

  void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                 NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;

  void
  SplitIntoTwoSpace(unsigned int k,
                    double space_best_score,
                    std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                    std::deque<std::shared_ptr<kTPMFQIterator>> &fq_stack,
                    std::deque<unsigned int> &iterate_times,
                    SpaceHeap *collection) override;

  // split fr which has only one node
  void SplitTopFq(unsigned int k,
                  double space_best_score,
                  std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                  std::deque<std::shared_ptr<kTPMFQIterator>> &fq_stack,
                  std::deque<unsigned int> &iterate_times,
                  SpaceHeap *collection);

  void
  MergeUpSubSpace(unsigned int k,
                  double space_best_score,
                  std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                  std::deque<std::shared_ptr<kTPMFQIterator>> &fq_stack,
                  std::deque<unsigned int> &iterate_times,
                  SpaceHeap *collection) override;
};

class kTPMOWIterator: public kTPMList{
 private:
  std::vector<kTPMElement> pool_;
 public:
  kTPMOWIterator() = default;
  ~kTPMOWIterator(){};
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
                    std::deque<std::shared_ptr<kTPMFQIterator>> &fq_stack,
                    std::deque<unsigned int> &iterate_times,
                    SpaceHeap *collection) override;

  void
  MergeUpSubSpace(unsigned int k,
                  double space_best_score,
                  std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                  std::deque<std::shared_ptr<kTPMFQIterator>> &fq_stack,
                  std::deque<unsigned int> &iterate_times,
                  SpaceHeap *collection)override;
};

// FQ may have cost itself
class kTPMFQIterator: public kTPMList{
 private:
  double node_score_;
  std::vector<std::shared_ptr<kTPMList>> fr_ow_iterators_;
  std::vector<NodeOneChildVarPredicatesPtr> types_predicates_;
  TYPE_ENTITY_LITERAL_ID node_id_;
 public:

  virtual ~kTPMFQIterator(){};
  void AddOneTypePredicate(NodeOneChildVarPredicatesPtr p){this->types_predicates_.push_back(p);}
  explicit kTPMFQIterator(int k, TYPE_ENTITY_LITERAL_ID node_id, int child_type_num, double node_score):
      node_score_(node_score),node_id_(node_id)
  {   this->fr_ow_iterators_.reserve(child_type_num);};
  void GetFirst(unsigned int k) override;
  inline TYPE_ENTITY_LITERAL_ID GetNodeID(){return this->node_id_;};

  void Insert(std::vector<std::shared_ptr<kTPMList>> FR_OW_iterators);
  void Insert(std::shared_ptr<kTPMList> FR_OW_iterator);
  inline unsigned int ChildTypeNum() {return this->fr_ow_iterators_.size();}
  virtual void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
  virtual OrderedListType Type() override {return OrderedListType::FQ;};
  // 谓词这里的ID还没搞好
  void
  SplitIntoTwoSpace(unsigned int k,
                    double space_best_score,
                    std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                    std::deque<std::shared_ptr<kTPMFQIterator>> &fq_stack,
                    std::deque<unsigned int> &iterate_times,
                    SpaceHeap *collection) override;
  void
  MergeUpSubSpace(unsigned int k,
                  double space_best_score,
                  std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                  std::deque<std::shared_ptr<kTPMFQIterator>> &fq_stack,
                  std::deque<unsigned int> &iterate_times,
                  SpaceHeap *collection)override;
};

class Subspace{
 public:
  std::vector<TYPE_ENTITY_LITERAL_ID> one_node_space_;
  std::deque<std::shared_ptr<kTPMFQIterator>> fq_stack_;
  // records which i-th child of fq_stack[i] is being visited
  std::deque<unsigned int> iterate_times_;
  // end_node_ must be FR/OW
  std::shared_ptr<kTPMList> end_node_;
  double best_score_;
  bool operator<(const Subspace& other) const
  {
    return this->best_score_ < other.best_score_;
  }
};

class SubspaceLess{
 public:
  inline bool operator()(const std::unique_ptr<Subspace>& p1,const std::unique_ptr<Subspace>& p2) const
  {
    return p1->best_score_<p2->best_score_;
  }
};


#endif //TOPK_SUBSPACE_H_
