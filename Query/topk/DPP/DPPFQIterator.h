
#ifndef PRIVATETOPK_QUERY_TOPK_DPP_DPPFQITERATOR_H_
#define PRIVATETOPK_QUERY_TOPK_DPP_DPPFQITERATOR_H_
#include "./DPPList.h"
#include "./GlobalQueue.h"

class DPPFRIterator;
// FQ may have cost itself
class DPPFQIterator: public DPPList{
 private:
  double node_score_;
  minmax::MinMaxHeap<FqElement> queue_;
  std::vector<std::shared_ptr<DPPList>> fr_ow_iterators_;
  ePool e_pool_;
  DynamicTrie dynamic_trie_;
  std::vector<NodeOneChildVarPredicatesPtr> types_predicates_;
  TYPE_ENTITY_LITERAL_ID node_id_;
  size_t child_type_num_;
  size_t branch_counter_;
  double cost1_;
  std::vector<double> sin_seq_cost_;
  std::vector< std::pair<DPPFRIterator*,unsigned int> > temp_parents;
  std::set< DPPFRIterator*> parents_;
  void GQ_Update(GlobalQueue *global_queue,std::unique_ptr<FqElement> e);
 public:
  bool visited_;
  void InsertEPool(FqElement e){
    this->pool_.push_back(element{this->node_id_,static_cast<unsigned int>(this->pool_.size()),e.cost});
    this->e_pool_.push_back(std::move(e));
  }
  void IncreaseBranchNum() { branch_counter_++;}
  std::size_t EPoolSize(){return this->e_pool_.size();}
  void TryTriggeringRootSeq(size_t j_type,GlobalQueue *global_queue);
  void TriggerASingletonSeq(unsigned int k,size_t j_type,GlobalQueue *global_queue);
  virtual ~DPPFQIterator(){};
  std::shared_ptr<DPPList> GetChild(unsigned int i_th){return this->fr_ow_iterators_[i_th];}
  void AddOneTypePredicate(NodeOneChildVarPredicatesPtr p){this->types_predicates_.push_back(p);}
  explicit DPPFQIterator(int k, TYPE_ENTITY_LITERAL_ID node_id, int child_type_num, double node_score,int fq_var):
      node_score_(node_score), dynamic_trie_(child_type_num,k),node_id_(node_id),branch_counter_(0),cost1_(0.0),
      child_type_num_(child_type_num),sin_seq_cost_(child_type_num,0.0),visited_(false)
  {
    this->var_id = fq_var;
    this->fr_ow_iterators_.reserve(child_type_num);
  };
  void GetFirst(unsigned int k,GlobalQueue *global_queue);
  void TryGetNext(unsigned int k,GlobalQueue *global_queue);
  TYPE_ENTITY_LITERAL_ID GetNodeID(){return this->node_id_;};

  void Insert(std::vector<std::shared_ptr<DPPList>> FR_OW_iterators);
  void Insert(std::shared_ptr<DPPList> FR_OW_iterator);
  bool NextEPoolElement(unsigned int k, std::shared_ptr<DPPList> node_pointer,
                               unsigned int index,unsigned int j_th_child,GlobalQueue *global_queue);
  bool NextIPoolElement(unsigned int k , size_t j_th_type, unsigned int index, GlobalQueue *global_queue);
  static double DeltaCost(std::shared_ptr<DPPList> FR_OW_iterator, int index);
  virtual void GetResult(int i_th,std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
  virtual OrderedListType Type() override {return OrderedListType::FQ;};
  std::vector< std::pair<DPPFRIterator*,unsigned int> >& TempParents(){return temp_parents;}
  std::set< DPPFRIterator*>& Parents(){return this->parents_;}
  void CleanTempParents(){ this->temp_parents.clear();}
  void AddParent(DPPFRIterator* fq_iterator,unsigned int index);
  TYPE_ENTITY_LITERAL_ID GetID() {return this->node_id_;}
};


#endif //PRIVATETOPK_QUERY_TOPK_DPP_DPPFQITERATOR_H_
