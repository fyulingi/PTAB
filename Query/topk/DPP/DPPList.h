
#ifndef PRIVATETOPK_QUERY_TOPK_DPP_PLIST_H_
#define PRIVATETOPK_QUERY_TOPK_DPP_PLIST_H_

#include "../../../Util/Util.h"
#include "../CommonDP/Pool.h"
#include "../CommonDP/MinMaxHeap.hpp"
#include "../CommonDP/DynamicTrie.h"

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
class DPPList{
 public:
  int var_id;
  Pool pool_;
  size_t GQ_pos_;
  DPPList():GQ_pos_(-1),var_id(-1){};
  // predicates_vec[i] is the predicates of child[i]
  // FQs and FRs will use this field.
  virtual OrderedListType Type(){return OrderedListType::UnDefined;};
  virtual void GetResult(int i_th,std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                         NodeOneChildVarPredicatesPtr predicate_information = nullptr)=0;
  virtual ~DPPList(){};
};


#endif //PRIVATETOPK_QUERY_TOPK_DPP_PLIST_H_
