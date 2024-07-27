
#include "Take2Subspace.h"

// getFirst() also implemented here
void Take2FRIterator::GetFirst(unsigned int k) {
  return;
}


/**
 * Insert One FQ iterator into the FR, and may change the
 * top Take2Element of the FR
 * @param k
 * @param fq_pointer
 */
void Take2FRIterator::Insert(unsigned int k,
                           TYPE_ENTITY_LITERAL_ID fq_id,
                           std::shared_ptr<Take2List> fq_pointer,
                           OnePointPredicatePtr predicates_vec) {
  this->heap_.emplace_back(fq_id,
                           fq_pointer->top_element_ptr_->cost,
                           std::move(fq_pointer),
                           std::move(predicates_vec));
}

void Take2FRIterator::MakeHeap()
{
  std::make_heap(this->heap_.begin(),this->heap_.end());
  this->top_element_ptr_ = std::unique_ptr<Take2frElement>(new Take2frElement(this->heap_[0].node,
                                                                              this->heap_[0].cost,
                                                                              this->heap_[0].fq_,
                                                                              this->heap_[0].predicate_));
  this->Take2List::top_element_ptr_ = std::unique_ptr<Take2Element>(new Take2Element(this->heap_[0].node,this->heap_[0].cost));
}

/**
 *
 * @param i_th
 * @param record
 * @param predicate_information not used, because FR itself saves the information
 */
void Take2FRIterator::GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                           NodeOneChildVarPredicatesPtr predicate_information) {
  auto &fq = *this->top_element_ptr_;
  auto fq_id = fq.node;
  auto fq_pointer = this->top_element_ptr_->fq_;
  // auto fq_pointer = this->fqs_map_[fq_id];
  // if(this->type_predicates_->find(fq_id)!=this->type_predicates_->end()) {
  if(this->top_element_ptr_->predicate_!= nullptr){
    for (auto pre_id:*this->top_element_ptr_->predicate_)
      record->push_back(pre_id);
  }
  fq_pointer->GetResult(record);
}

void
Take2FRIterator::SplitIntoThreeSpace(unsigned int k,
                                     double space_best_score,
                                     std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                                     std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                                     std::deque<unsigned int> &iterate_times,
                                     Take2SpaceHeap *collection) {
  auto top_id = this->top_element_ptr_->node;
  auto top_cost_ = this->top_element_ptr_->cost;
  auto remaining_cost = space_best_score - top_cost_;

  auto top_fr = std::make_shared<Take2FRIterator>();


  auto top_fq = this->top_element_ptr_->fq_;
  auto top_predicate =this->top_element_ptr_->predicate_;

  top_fr->Insert(k,top_id,top_fq,top_predicate);
  top_fr->MakeHeap();

  auto nums = this->heap_.size();
  bool first_child_valid = nums >1;
  bool second_child_valid = nums >2;

  auto first_fr = std::make_shared<Take2FRIterator>();
  auto second_fr = std::make_shared<Take2FRIterator>();

  unsigned int p = 1;
  unsigned int step = 1;
  while(p<nums)
  {
    // add to the left fr
    for(int j=0;j<step && p+j<nums;j++)
      first_fr->Insert(this->heap_[p+j]);
    p += step;
    // add to the right fr
    for(int j=0;j<step && p+j<nums;j++)
      second_fr->Insert(this->heap_[p+j]);
    p += step;
    step = step + step;
  }

  if(first_child_valid)
  {
    first_fr->top_element_ptr_ = std::unique_ptr<Take2frElement>(new Take2frElement(first_fr->heap_[0].node,
                                                                                    first_fr->heap_[0].cost,
                                                                                    first_fr->heap_[0].fq_,
                                                                                    first_fr->heap_[0].predicate_));
    first_fr->Take2List::top_element_ptr_ = std::unique_ptr<Take2Element>(new Take2Element(this->heap_[0].node,this->heap_[0].cost));
    auto new_sub_space = std::unique_ptr<Take2Subspace>(new Take2Subspace);
    new_sub_space->fq_stack_ = fq_stack;
    new_sub_space->iterate_times_ = iterate_times;
    new_sub_space->end_node_ = first_fr;
    new_sub_space->one_node_space_ = record;
    new_sub_space->best_score_ = remaining_cost + first_fr->top_element_ptr_->cost;
    collection->push(std::move(new_sub_space));
  }

  if(second_child_valid)
  {
    second_fr->top_element_ptr_ = std::unique_ptr<Take2frElement>(new Take2frElement(second_fr->heap_[0].node,
                                                                                     second_fr->heap_[0].cost,
                                                                                     second_fr->heap_[0].fq_,
                                                                                     second_fr->heap_[0].predicate_));
    second_fr->Take2List::top_element_ptr_ = std::unique_ptr<Take2Element>(new Take2Element(this->heap_[0].node,this->heap_[0].cost));
    auto new_sub_space = std::unique_ptr<Take2Subspace>(new Take2Subspace);
    new_sub_space->fq_stack_ = fq_stack;
    new_sub_space->iterate_times_ = iterate_times;
    new_sub_space->end_node_ = second_fr;
    new_sub_space->one_node_space_ = record;
    new_sub_space->best_score_ = remaining_cost + second_fr->top_element_ptr_->cost;
    collection->push(std::move(new_sub_space));
  }

  top_fr->SplitTopFq(k,space_best_score,
                     record,fq_stack,iterate_times,collection);
}

Take2FRIterator::Take2FRIterator() {
  // this->type_predicates_=std::make_shared<NodeOneChildVarPredicates>();
}

void Take2FRIterator::SplitTopFq(unsigned int k,
                                double space_best_score,
                                std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                                std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                                std::deque<unsigned int> &iterate_times,
                                Take2SpaceHeap *collection){
  auto fq_id = this->top_element_ptr_->node;
  // auto fq = this->fqs_map_[fq_id];
  auto fq = this->top_element_ptr_->fq_;
  if(this->top_element_ptr_->predicate_!= nullptr){
    for (auto pre_id:*this->top_element_ptr_->predicate_)
      record.push_back(pre_id);
  }
//  if(this->type_predicates_->find(fq_id)!=this->type_predicates_->end()) {
//    auto predicates = (*this->type_predicates_)[fq_id];
//    for (auto pre_id:*predicates)
//      record.push_back(pre_id);
//  }
  fq->SplitIntoThreeSpace(k, space_best_score, record, fq_stack, iterate_times, collection);
}

void Take2FRIterator::MergeUpSubSpace(unsigned int k,
                                     double space_best_score,
                                     std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                                     std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                                     std::deque<unsigned int> &iterate_times,
                                     Take2SpaceHeap *collection) {
  this->SplitIntoThreeSpace(k, space_best_score, record, fq_stack, iterate_times, collection);
  if(!fq_stack.empty()) {
    auto &back_fq = fq_stack.back();
    back_fq->MergeUpSubSpace(k, space_best_score, record, fq_stack, iterate_times, collection);
  }
}

void Take2FRIterator::Insert(Take2frElement &ele) {
  this->heap_.push_back(ele);
}

void Take2OWIterator::GetFirst(unsigned int k) {
}

void Take2OWIterator::Insert(unsigned int k,
    const std::vector<TYPE_ENTITY_LITERAL_ID>& ids,
                           const std::vector<double>& scores) {
//  struct ScorePair{
//    TYPE_ENTITY_LITERAL_ID id;
//    double cost;
//    bool operator<(const ScorePair& other) const{return this->cost<other.cost;};
//  };
//  std::vector<ScorePair> ranks;
//  ranks.reserve(ids.size());
//  for(unsigned int  i=0;i<ids.size();i++)
//    ranks.push_back(ScorePair{ids[i],scores[i]});
//  std::sort(ranks.begin(),ranks.end());
//  for(unsigned int i=0;i<k && i < ranks.size();i++)
//  {
//    Take2Element e{};
//    e.cost = ranks[i].cost;
//    e.node = ranks[i].id;
//    this->pool_.push_back(e);
//  }
  for (unsigned int i = 0; i < ids.size(); i++)
    this->heap_.emplace_back(ids[i],scores[i]);
}

/**
 * Insert a made heap
 * @param k
 * @param ids
 * @param scores
 */
void Take2OWIterator::RankedInsert(unsigned int k,
                            const std::vector<TYPE_ENTITY_LITERAL_ID>& ids,
                            const std::vector<double>& scores)
{
  for (unsigned int i = 0; i < ids.size(); i++)
    this->heap_.emplace_back(ids[i],scores[i]);
}

/**
 * Default score : 0.0
 * @param ids
 */
void Take2OWIterator::Insert(unsigned int k,
                           const std::vector<TYPE_ENTITY_LITERAL_ID> &ids) {
  for (unsigned int i = 0; i < ids.size(); i++)
    this->heap_.emplace_back(ids[i],0.0);
}

void Take2OWIterator::Insert(Take2Element &ele)
{
  this->heap_.push_back(ele);
}
void Take2OWIterator::MakeHeap()
{
  std::make_heap(this->heap_.begin(),this->heap_.end());
  this->top_element_ptr_ = std::unique_ptr<Take2Element>(new Take2Element(this->heap_[0]));
}
void Take2OWIterator::GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                              NodeOneChildVarPredicatesPtr predicate_information) {
  auto &i_th_Take2Element = *this->top_element_ptr_;
  auto node_id = i_th_Take2Element.node;
  auto predicates = (*predicate_information)[node_id];
  for(auto predicate_id:*predicates)
    record->push_back(predicate_id);
  record->push_back(node_id);
}

void
Take2OWIterator::SplitIntoThreeSpace(unsigned int k,
                                     double space_best_score,
                                     std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                                     std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                                     std::deque<unsigned int> &iterate_times,
                                     Take2SpaceHeap *collection) {
  auto top_cost_ = this->heap_[0].cost;
  auto remaining_cost = space_best_score - top_cost_;

  auto nums = this->heap_.size();
  bool first_child_valid = nums >1;
  bool second_child_valid = nums >2;

  auto first_ow = std::make_shared<Take2OWIterator>();
  auto second_ow = std::make_shared<Take2OWIterator>();

  unsigned int p = 1;
  unsigned int step = 1;
  while(p<nums)
  {
    // add to the left fr
    for(int j=0;j<step && p+j<nums;j++)
      first_ow->Insert(this->heap_[p+j]);
    p += step;
    // add to the right fr
    for(int j=0;j<step && p+j<nums;j++)
      second_ow->Insert(this->heap_[p+j]);
    p += step;
    step = step + step;
  }

  if(first_child_valid)
  {
    first_ow->top_element_ptr_ = std::unique_ptr<Take2Element>(new Take2Element(first_ow->heap_[0]));
    auto new_sub_space = std::unique_ptr<Take2Subspace>(new Take2Subspace);
    new_sub_space->fq_stack_ = fq_stack;
    new_sub_space->iterate_times_ = iterate_times;
    new_sub_space->end_node_ = first_ow;
    new_sub_space->one_node_space_ = record;
    new_sub_space->best_score_ = remaining_cost + first_ow->top_element_ptr_->cost;
    collection->push(std::move(new_sub_space));
  }

  if(second_child_valid)
  {
    second_ow->top_element_ptr_ = std::unique_ptr<Take2Element>(new Take2Element(second_ow->heap_[0]));
    auto new_sub_space = std::unique_ptr<Take2Subspace>(new Take2Subspace);
    new_sub_space->fq_stack_ = fq_stack;
    new_sub_space->iterate_times_ = iterate_times;
    new_sub_space->end_node_ = second_ow;
    new_sub_space->one_node_space_ = record;
    new_sub_space->best_score_ = remaining_cost + second_ow->top_element_ptr_->cost;
    collection->push(std::move(new_sub_space));
  }
  record.push_back(this->heap_[0].node);
}

void
Take2OWIterator::MergeUpSubSpace(unsigned int k,
                                double space_best_score,
                                std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                                std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                                std::deque<unsigned int> &iterate_times,
                                Take2SpaceHeap *collection)
{
  this->SplitIntoThreeSpace(k, space_best_score, record, fq_stack, iterate_times, collection);
  if(!fq_stack.empty()) {
    auto &back_fq = fq_stack.back();
    back_fq->MergeUpSubSpace(k, space_best_score, record, fq_stack, iterate_times, collection);
  }
}

void Take2FQIterator::GetFirst(unsigned int k) {
  double cost = 0;
  for (unsigned int j = 0; j < this->fr_ow_iterators_.size(); j++) {
    cost += this->fr_ow_iterators_[j]->top_element_ptr_->cost;
  }
  auto ele_ptr = std::unique_ptr<Take2Element>(new Take2Element(this->node_id_,cost));
  this->top_element_ptr_ = std::move(ele_ptr);
}

void Take2FQIterator::Insert(std::shared_ptr<Take2List> FR_OW_iterator) {
  this->fr_ow_iterators_.push_back(std::move(FR_OW_iterator));
}

/**
 * Insert a bulk of FR or OW iterators.
 * inserting one certain type each time each time
 * certain type [i] specified by it's the i-th child of its father
 * @param FR_OW_iterators
 */
void Take2FQIterator::Insert(std::vector<std::shared_ptr<Take2List>> FR_OW_iterators) {
  this->fr_ow_iterators_ = std::move(FR_OW_iterators);
}

void Take2FQIterator::GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                               NodeOneChildVarPredicatesPtr predicate_information) {
  record->push_back(this->node_id_);
  for(unsigned int i =0; i<this->fr_ow_iterators_.size(); i++) {
    if(fr_ow_iterators_[i]->Type() ==OrderedListType::OW) {
      auto ow_predicates = this->types_predicates_[i];
      fr_ow_iterators_[i]->GetResult( record,ow_predicates);
    }
    else
      fr_ow_iterators_[i]->GetResult( record);
  }
}

void
Take2FQIterator::SplitIntoThreeSpace(unsigned int k,
                                     double space_best_score,
                                     std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                                     std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                                     std::deque<unsigned int> &iterate_times,
                                     Take2SpaceHeap *collection) {
  record.push_back(this->node_id_);
  fq_stack.push_back(this->shared_from_base<Take2FQIterator>());
  iterate_times.push_back(0);
  for(unsigned int i =0 ;i<this->fr_ow_iterators_.size();i++)
  {
    auto &fr_ow_iterator = this->fr_ow_iterators_[i];
    if(fr_ow_iterator->Type()==OrderedListType::OW) {
      auto ow_predicates = this->types_predicates_[i];
      auto node_id = fr_ow_iterator->top_element_ptr_->node;
      auto predicates = (*ow_predicates)[node_id];
      for (auto predicate_id:*predicates)
        record.push_back(predicate_id);
    }
    fr_ow_iterator->SplitIntoThreeSpace(k, space_best_score, record, fq_stack, iterate_times, collection);
    iterate_times.back()++;
  }
  iterate_times.pop_back();
  fq_stack.pop_back();
}

void Take2FQIterator::MergeUpSubSpace(unsigned int k,
                                    double space_best_score,
                                    std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                                    std::deque<std::shared_ptr<Take2FQIterator>> &fq_stack,
                                    std::deque<unsigned int> &iterate_times,
                                    Take2SpaceHeap *collection) {

  unsigned int times = iterate_times.back();
  for(unsigned int i=times+1;i<this->fr_ow_iterators_.size();i++)
  {
    auto &fr_ow_iterator = this->fr_ow_iterators_[i];
    if(fr_ow_iterator->Type()==OrderedListType::OW) {
      auto ow_predicates = this->types_predicates_[i];
      auto node_id = fr_ow_iterator->top_element_ptr_->node;
      auto predicates = (*ow_predicates)[node_id];
      for (auto predicate_id:*predicates)
        record.push_back(predicate_id);
    }
    iterate_times.back()++;
    this->fr_ow_iterators_[i]->SplitIntoThreeSpace(k, space_best_score, record, fq_stack, iterate_times, collection);
  }
  fq_stack.pop_back();
  iterate_times.pop_back();
  if(!fq_stack.empty())
  {
    auto upper_fq = fq_stack.back();
    upper_fq->MergeUpSubSpace(k,space_best_score,record,fq_stack,iterate_times,collection);
  }
}
