
#ifndef PRIVATETOPK_QUERY_TOPK_DPP_DPPOWITERATOR_H_
#define PRIVATETOPK_QUERY_TOPK_DPP_DPPOWITERATOR_H_

#include "./DPPList.h"
#include "./GlobalQueue.h"

class DPPOWIterator: public DPPList{
 private:
  size_t valid_num = 0;
  // sorting right after building, so don't need a heap to get top-k
  // the top-k result already in the pool
  // record a valid num, for DP-P
 public:
  DPPOWIterator(): valid_num(0){};
  virtual ~DPPOWIterator(){};
  OrderedListType Type() override {return OrderedListType::OW;};
  void TryGetNext(unsigned int k);
  size_t Size(){return  valid_num;};
  // Insert a bulk of gStore Node ids and their scores
  void Insert(unsigned int k, const std::vector<TYPE_ENTITY_LITERAL_ID>& ids, const std::vector<double>& scores);
  void Insert(unsigned int k,const std::vector<TYPE_ENTITY_LITERAL_ID> &ids);
  virtual void GetResult(int i_th,std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
};

#endif //PRIVATETOPK_QUERY_TOPK_DPP_DPPOWITERATOR_H_
