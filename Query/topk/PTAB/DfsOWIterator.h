
#ifndef PRIVATETOPK_QUERY_TOPK_Dfs_DfsOWITERATOR_H_
#define PRIVATETOPK_QUERY_TOPK_Dfs_DfsOWITERATOR_H_

#include "./DfsList.h"
#include "./HopIndex.h"
class DfsOWIterator: public DfsList{
 private:
  // sorting right after building, so don't need a heap to get top-k
  // the top-k result already in the pool
 public:
  DfsOWIterator()=default;
  virtual ~DfsOWIterator(){};
  OrderedListType Type() override {return OrderedListType::OW;};
  bool Empty() {
    return this->pool_.empty();
  };
  void TryGetNext(DfsHelpInfo *dfs_help_info);
  // Insert a bulk of gStore Node ids and their scores
  void Insert(unsigned int k, const std::vector<TYPE_ENTITY_LITERAL_ID>& ids, const std::vector<double>& scores);
  void Insert(unsigned int k,const std::vector<TYPE_ENTITY_LITERAL_ID> &ids);
  virtual void GetResult(int i_th,std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr) override;
  void Estimate() override { return; };
};

#endif //PRIVATETOPK_QUERY_TOPK_Dfs_DfsOWITERATOR_H_
