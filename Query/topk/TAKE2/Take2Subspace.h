
#ifndef TAKE2_SUBSPACE_H_
#define TAKE2_SUBSPACE_H_
#include "../../../Util/Util.h"
#include "../CommonDP/Pool.h"
#include "../CommonDP/MinMaxHeap.hpp"

class Take2Subspace;
class Take2SubspaceLess;
class Take2FQIterator;
using Take2SpaceHeap = minmax::MinMaxHeap<std::unique_ptr<Take2Subspace>,std::vector<std::unique_ptr<Take2Subspace>>,Take2SubspaceLess>;

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
  selecting top-k min Take2Element */
class Take2List: public std::enable_shared_from_this<Take2List>{
 public:
  std::unique_ptr<Take2Element> top_element_ptr_;
  Take2List()=default;
  // predicates_vec[i] is the predicates of child[i]
  // FQs and FRs will use this field.
  virtual OrderedListType Type(){return OrderedListType::UnDefined;};
  virtual void GetFirst(unsigned int k)=0;
  virtual void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr)=0;
  virtual ~Take2List()=default;

  // cut a node from this iterator and this Take2List will changed to two new List
  // the origin will remain
  virtual void
  SplitIntoThreeSpace(unsigned int k,
                      double space_best_score,
                      std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                      std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                      std::deque<unsigned int> &iterate_times,
                      Take2SpaceHeap *collection)=0;
  virtual void
  MergeUpSubSpace(unsigned int k,
                    double space_best_score,
                    std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                    std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                    std::deque<unsigned int> &iterate_times,
                    Take2SpaceHeap *collection)=0;

 protected:
  template <typename Derived>
  std::shared_ptr<Derived> shared_from_base()
  {
    return std::static_pointer_cast<Derived>(shared_from_this());
  }
};

struct Take2frElement{
  TYPE_ENTITY_LITERAL_ID node;
  double cost;
  std::shared_ptr<Take2List> fq_;
  OnePointPredicatePtr predicate_;
  Take2frElement(TYPE_ENTITY_LITERAL_ID node_id,double subtree_cost,
                 std::shared_ptr<Take2List> fq,OnePointPredicatePtr predicate):
      node(node_id),cost(subtree_cost),fq_(std::move(fq)),predicate_(std::move(predicate))
  {}
  // for max heap
  bool operator < (const Take2frElement &other) const {return this->cost > other.cost;}
  bool operator!=(const Take2frElement &other) const {
    return this->node!=other.node || this->cost!=other.cost;
  }
};

class Take2FRIterator: public Take2List {
 private:
  std::vector<Take2frElement> heap_;
  // std::vector<TYPE_ENTITY_LITERAL_ID> fq_ids_;
  // std::vector<OnePointPredicatePtr> predicates_;
  // std::vector<std::shared_ptr<Take2List>> fqs_;

 public:
  std::unique_ptr<Take2frElement> top_element_ptr_;
  Take2FRIterator();
  ~Take2FRIterator(){};
  OrderedListType Type() override {return OrderedListType::FR;};
  void GetFirst(unsigned int k) override;

  // Insert a bulk of FQ iterator, all the same type
  void Insert(unsigned int k,
              TYPE_ENTITY_LITERAL_ID fq_id,
              std::shared_ptr<Take2List> FQ_iterator,
              OnePointPredicatePtr predicates_vec);
  void Insert(Take2frElement &ele);
  void MakeHeap();
  void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                 NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;

  void
  SplitIntoThreeSpace(unsigned int k,
                      double space_best_score,
                      std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                      std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                      std::deque<unsigned int> &iterate_times,
                      Take2SpaceHeap *collection) override;

  // split fr which has only one node
  void SplitTopFq(unsigned int k,
                  double space_best_score,
                  std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                  std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                  std::deque<unsigned int> &iterate_times,
                  Take2SpaceHeap *collection);

  void
  MergeUpSubSpace(unsigned int k,
                  double space_best_score,
                  std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                  std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                  std::deque<unsigned int> &iterate_times,
                  Take2SpaceHeap *collection) override;
};

class Take2OWIterator: public Take2List{
 private:
  std::vector<Take2Element> heap_;
 public:
  Take2OWIterator() = default;
  ~Take2OWIterator(){};
  OrderedListType Type() override {return OrderedListType::OW;};
  void GetFirst(unsigned int k) override;

  // Insert a bulk of gStore Node ids and their scores
  void Insert(unsigned int k, const std::vector<TYPE_ENTITY_LITERAL_ID>& ids, const std::vector<double>& scores);
  void Insert(unsigned int k,const std::vector<TYPE_ENTITY_LITERAL_ID> &ids);
  void Insert(Take2Element &ele);
  void RankedInsert(unsigned int k, const std::vector<TYPE_ENTITY_LITERAL_ID>& ids, const std::vector<double>& scores);
  void MakeHeap();
  void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
  void
  SplitIntoThreeSpace(unsigned int k,
                      double space_best_score,
                      std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                      std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                      std::deque<unsigned int> &iterate_times,
                      Take2SpaceHeap *collection) override;

  void
  MergeUpSubSpace(unsigned int k,
                  double space_best_score,
                  std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                  std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                  std::deque<unsigned int> &iterate_times,
                  Take2SpaceHeap *collection)override;
};

// FQ may have cost itself
class Take2FQIterator: public Take2List{
 private:
  double node_score_;
  std::vector<std::shared_ptr<Take2List>> fr_ow_iterators_;
  std::vector<NodeOneChildVarPredicatesPtr> types_predicates_;
  TYPE_ENTITY_LITERAL_ID node_id_;
 public:

  virtual ~Take2FQIterator(){};
  void AddOneTypePredicate(NodeOneChildVarPredicatesPtr p){this->types_predicates_.push_back(p);}
  explicit Take2FQIterator(int k, TYPE_ENTITY_LITERAL_ID node_id, int child_type_num, double node_score):
      node_score_(node_score),node_id_(node_id)
  {   this->fr_ow_iterators_.reserve(child_type_num);};
  void GetFirst(unsigned int k) override;
  inline TYPE_ENTITY_LITERAL_ID GetNodeID(){return this->node_id_;};

  void Insert(std::vector<std::shared_ptr<Take2List>> FR_OW_iterators);
  void Insert(std::shared_ptr<Take2List> FR_OW_iterator);
  inline unsigned int ChildTypeNum() {return this->fr_ow_iterators_.size();}
  virtual void GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
  virtual OrderedListType Type() override {return OrderedListType::FQ;};
  void
  SplitIntoThreeSpace(unsigned int k,
                      double space_best_score,
                      std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                      std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                      std::deque<unsigned int> &iterate_times,
                      Take2SpaceHeap *collection) override;
  void
  MergeUpSubSpace(unsigned int k,
                  double space_best_score,
                  std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                  std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                  std::deque<unsigned int> &iterate_times,
                  Take2SpaceHeap *collection)override;
};

class Take2Subspace{
 public:
  std::vector<TYPE_ENTITY_LITERAL_ID> one_node_space_;
  std::deque<std::shared_ptr<Take2FQIterator>> fq_stack_;
  // records which i-th child of fq_stack[i] is being visited
  std::deque<unsigned int> iterate_times_;
  // end_node_ must be FR/OW
  std::shared_ptr<Take2List> end_node_;
  double best_score_;
  bool operator<(const Take2Subspace& other) const
  {
    return this->best_score_ < other.best_score_;
  }
};

class Take2SubspaceLess{
 public:
  inline bool operator()(const std::unique_ptr<Take2Subspace>& p1,const std::unique_ptr<Take2Subspace>& p2) const
  {
    return p1->best_score_<p2->best_score_;
  }
};


#endif //TAKE2_SUBSPACE_H_
