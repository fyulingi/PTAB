
#ifndef PRIVATETOPK_QUERY_TOPK_DPP_DPPFRITERATOR_H_
#define PRIVATETOPK_QUERY_TOPK_DPP_DPPFRITERATOR_H_

#include "./DPPList.h"
#include "./GlobalQueue.h"

class DPPFRIterator: public DPPList {
 private:
  NodeOneChildVarPredicatesPtr type_predicates_;
  minmax::MinMaxHeap<element> queue_;
  std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DPPFQIterator>> fqs_map_;
  void GQ_Update(GlobalQueue *global_queue,element e, unsigned i_th_type);
  DPPFQIterator *parent_fq_;
 public:
  DPPFRIterator();
  virtual ~DPPFRIterator(){};
  OrderedListType Type() override {return OrderedListType::FR;};

  void GetFirst();

  void TryGetNext(unsigned int k,GlobalQueue *global_queue,
                  DPPFQIterator* fq_iterator,unsigned int i_th_child);

  // Insert a bulk of FQ iterator, all the same type
  void Insert(TYPE_ENTITY_LITERAL_ID fq_id,
              std::shared_ptr<DPPFQIterator> fq_pointer,
              OnePointPredicatePtr predicates_vec);

  void Register(TYPE_ENTITY_LITERAL_ID fq_id,
              std::shared_ptr<DPPFQIterator> fq_pointer,
              OnePointPredicatePtr predicates_vec);

  void dInsert(unsigned int k,
               DPPFQIterator* child_fq_pointer,
               unsigned i_th_type,
               GlobalQueue *global_queue);

  void rInsert(unsigned int k,
               int index,
               DPPList* fq_pointer,
               unsigned i_th_type,
               GlobalQueue *global_queue);

  static double DeltaCost(std::shared_ptr<DPPFQIterator> node_pointer, int index);
  bool NextEPoolElement(unsigned int k,
                               std::shared_ptr<DPPFQIterator> node_pointer,
                               unsigned int index,
                               GlobalQueue *global_queue);
  void GetResult(int i_th,std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
  void SetParent(DPPFQIterator* parent){ this->parent_fq_ = parent; };

  std::map<TYPE_ENTITY_LITERAL_ID,std::shared_ptr<DPPFQIterator>>& GetFqMap(){ return this->fqs_map_;}
};



#endif //PRIVATETOPK_QUERY_TOPK_DPP_DPPFRITERATOR_H_
