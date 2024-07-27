
#include "SAESpaceSAE.h"

// getFirst() also implemented here
void SAEFRIterator::GetFirst(unsigned int k) {
  return;
}


/**
 * Insert One FQ iterator into the FR, and may change the
 * top SAEElement of the FR
 * @param k
 * @param fq_pointer
 */
void SAEFRIterator::Insert(unsigned int k,
                           TYPE_ENTITY_LITERAL_ID fq_id,
                           std::shared_ptr<SAEList> fq_pointer,
                           OnePointPredicatePtr predicates_vec) {
  //(*this->type_predicates_)[fq_id] = predicates_vec;
  //this->fqs_map_[fq_id] = fq_pointer;

  bool replace = false;
  if(this->top_element_ptr_ == nullptr)
    replace = true;
  else{
    auto cost = fq_pointer->top_element_ptr_->cost;
    if(cost<this->top_element_ptr_->cost)
      replace=true;
  }
  if(replace)
  {
    if(this->top_element_ptr_ == nullptr) {
      // write to the father class
      this->SAEList::top_element_ptr_ = std::unique_ptr<SAEElement>(new SAEElement);
      this->SAEList::top_element_ptr_->cost =  fq_pointer->top_element_ptr_->cost;
      this->SAEList::top_element_ptr_->node =  fq_id;
      // write my self
      auto e_ptr = std::unique_ptr<SAEfrElement>(new SAEfrElement);
      e_ptr->cost = fq_pointer->top_element_ptr_->cost;
      e_ptr->node = fq_id;
      e_ptr->fq_= std::move(fq_pointer);
      e_ptr->predicate_ = std::move(predicates_vec);
      this->top_element_ptr_ = std::move(e_ptr);

    }
    else
    {
      // save the old
      this->fq_ids_.push_back(top_element_ptr_->node);
      this->fqs_.push_back(top_element_ptr_->fq_);
      this->predicates_.push_back(top_element_ptr_->predicate_);
      // write the new
      this->SAEList::top_element_ptr_->cost =  fq_pointer->top_element_ptr_->cost;
      this->SAEList::top_element_ptr_->node =  fq_id;
      // write to this field
      this->top_element_ptr_->cost = fq_pointer->top_element_ptr_->cost;
      this->top_element_ptr_->node = fq_id;
      this->top_element_ptr_->fq_ = std::move(fq_pointer);
      this->top_element_ptr_->predicate_ = std::move(predicates_vec);
    }
  }
  else
  {
    this->fq_ids_.push_back(fq_id);
    this->fqs_.push_back(std::move(fq_pointer));
    this->predicates_.push_back(std::move(predicates_vec));
  }
}

/**
 *
 * @param i_th
 * @param record
 * @param predicate_information not used, because FR itself saves the information
 */
void SAEFRIterator::GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
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
SAEFRIterator::SplitIntoTwoSpace(unsigned int k,
                                  double space_best_score,
                                  std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                                  std::deque<std::shared_ptr<SAEFQIterator>> &fq_stack,
                                  std::deque<unsigned int> &iterate_times,
                                  SpaceHeap *collection) {
  auto top_id = this->top_element_ptr_->node;
  auto top_cost_ = this->top_element_ptr_->cost;
  auto remaining_cost = space_best_score - top_cost_;

  auto top_fr = std::make_shared<SAEFRIterator>();
  auto the_other_fr = std::make_shared<SAEFRIterator>();

  auto top_fq = this->top_element_ptr_->fq_;
  auto top_predicate =this->top_element_ptr_->predicate_;

  top_fr->Insert(k,top_id,top_fq,top_predicate);

//  for(unsigned int i=0;i<this->fq_ids_.size();i++)
//  {
//    auto fq_id = fq_pair.first;
//    if(fq_id == top_id)
//      continue;
//    auto fq_predicate = (*this->type_predicates_)[fq_id];
//    the_other_fr->Insert(k,fq_id,fq_pair.second,fq_predicate);
//  }
  bool the_other_empty = this->fq_ids_.empty();
  for(unsigned int i=0;i<this->fq_ids_.size();i++)
    the_other_fr->Insert(k, this->fq_ids_[i],
                         this->fqs_[i],
                         this->predicates_[i]);

  if(!the_other_empty)
  {
    auto new_sub_space = std::unique_ptr<SpaceSAE>(new SpaceSAE);
    new_sub_space->fq_stack_ = fq_stack;
    new_sub_space->iterate_times_ = iterate_times;
    new_sub_space->end_node_ = the_other_fr;
    new_sub_space->one_node_space_ = record;
    new_sub_space->best_score_ = remaining_cost + the_other_fr->top_element_ptr_->cost;
    collection->push(std::move(new_sub_space));
  }

  //if(!the_other_empty)
  //  the_other_fr->SplitTopFq(k,remaining_cost + the_other_fr->top_element_ptr_->cost,
  //                           record,fq_stack,iterate_times,collection);

  top_fr->SplitTopFq(k,space_best_score,
                     record,fq_stack,iterate_times,collection);
}

SAEFRIterator::SAEFRIterator() {
  // this->type_predicates_=std::make_shared<NodeOneChildVarPredicates>();
}

void SAEFRIterator::SplitTopFq(unsigned int k,
                                double space_best_score,
                                std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                                std::deque<std::shared_ptr<SAEFQIterator>> &fq_stack,
                                std::deque<unsigned int> &iterate_times,
                                SpaceHeap *collection){
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
  fq->SplitIntoTwoSpace(k,space_best_score,record,fq_stack,iterate_times,collection);
}

void SAEFRIterator::MergeUpSpaceSAE(unsigned int k,
                                     double space_best_score,
                                     std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                                     std::deque<std::shared_ptr<SAEFQIterator>> &fq_stack,
                                     std::deque<unsigned int> &iterate_times,
                                     SpaceHeap *collection) {
  this->SplitIntoTwoSpace(k,space_best_score,record,fq_stack,iterate_times,collection);
  if(!fq_stack.empty()) {
    auto &back_fq = fq_stack.back();
    back_fq->MergeUpSpaceSAE(k, space_best_score, record, fq_stack, iterate_times, collection);
  }
}

void SAEOWIterator::GetFirst(unsigned int k) {
}

void SAEOWIterator::Insert(unsigned int k,
    const std::vector<TYPE_ENTITY_LITERAL_ID>& ids,
                           const std::vector<double>& scores)
{
  struct ScorePair{
    TYPE_ENTITY_LITERAL_ID id;
    double cost;
    bool operator<(const ScorePair& other) const{return this->cost<other.cost;};
  };
  std::vector<ScorePair> ranks;
  this->pool_.reserve(k);
  ranks.reserve(ids.size());
  for(unsigned int  i=0;i<ids.size();i++)
    ranks.push_back(ScorePair{ids[i],scores[i]});
  std::sort(ranks.begin(),ranks.end());
  for(unsigned int i=0;i<k && i < ranks.size();i++)
  {
    SAEElement e{};
    e.cost = ranks[i].cost;
    e.node = ranks[i].id;
    this->pool_.push_back(e);
  }
  this->top_element_ptr_ = std::unique_ptr<SAEElement>(new SAEElement(this->pool_[0]));
}

void SAEOWIterator::RankedInsert(unsigned int k,
                            const std::vector<TYPE_ENTITY_LITERAL_ID>& ids,
                            const std::vector<double>& scores)
{
  for(unsigned int i=0;i<k && i < ids.size();i++)
  {
    SAEElement e{};
    e.cost = scores[i];
    e.node = ids[i];
    this->pool_.push_back(e);
  }
  this->top_element_ptr_ = std::unique_ptr<SAEElement>(new SAEElement(this->pool_[0]));
}

void SAEOWIterator::RankedInsert(unsigned int k,
                                  TYPE_ENTITY_LITERAL_ID id,
                                  double score)
{
  TakeAllElement e{};
  e.cost = score;
  e.node = id;
  this->pool_.push_back(e);
  if(this->pool_.size()==1)
    this->top_element_ptr_ = std::unique_ptr<TakeAllElement>(new TakeAllElement(this->pool_[0]));
}

/**
 * Default score : 0.0
 * @param ids
 */
void SAEOWIterator::Insert(unsigned int k,
                           const std::vector<TYPE_ENTITY_LITERAL_ID> &ids) {
  pool_.reserve(k);
  for(unsigned int i=0;i<k && i < ids.size();i++)
  {
    SAEElement e{};
    e.cost = 0.0;
    e.node = ids[i];
    this->pool_.push_back(e);
  }
  this->top_element_ptr_ = std::unique_ptr<SAEElement>(new SAEElement(this->pool_[0]));
}

void SAEOWIterator::GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                              NodeOneChildVarPredicatesPtr predicate_information) {
  auto &i_th_SAEElement = *this->top_element_ptr_;
  auto node_id = i_th_SAEElement.node;
  auto predicates = (*predicate_information)[node_id];
  for(auto predicate_id:*predicates)
    record->push_back(predicate_id);
  record->push_back(node_id);
}

void
SAEOWIterator::SplitIntoTwoSpace(unsigned int k,
                                  double space_best_score,
                                  std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                                  std::deque<std::shared_ptr<SAEFQIterator>> &fq_stack,
                                  std::deque<unsigned int> &iterate_times,
                                  SpaceHeap *collection) {
  auto top_cost_ = this->pool_[0].cost;
  auto remaining_cost = space_best_score - top_cost_;

  auto the_other_ow = std::make_shared<SAEOWIterator>();

//  bool no_coefficient = std::all_of(this->pool_.cbegin(), this->pool_.cend(),
//                                    [](const SAEElement &ele){ return ele.cost== 0.0; });
//
//
//  std::vector<TYPE_ENTITY_LITERAL_ID> ids;
//  std::vector<double> scores;
//  ids.reserve(this->pool_.size()-1);
//  scores.reserve(this->pool_.size()-1);
//
//  for(int i=1;i<this->pool_.size();i++)
//  {
//    ids.push_back(this->pool_[i].node);
//    scores.push_back(this->pool_[i].cost);
//  }

  bool the_other_empty = this->pool_.size() == 1;
  if(!the_other_empty) {
//    if (no_coefficient)
//      the_other_ow->Insert(k, ids);
//    else
    for(unsigned int i = 1;i < this->pool_.size();i++)
      the_other_ow->RankedInsert(k, this->pool_[i].node, this->pool_[i].cost);
    auto new_sub_space = std::unique_ptr<SpaceSAE>(new SpaceSAE);
    new_sub_space->fq_stack_ = fq_stack;
    new_sub_space->iterate_times_ = iterate_times;
    new_sub_space->end_node_ = the_other_ow;
    new_sub_space->one_node_space_ = record;
    new_sub_space->best_score_ = remaining_cost + the_other_ow->top_element_ptr_->cost;
    collection->push(std::move(new_sub_space));
  }
  record.push_back(this->pool_[0].node);
}

void
SAEOWIterator::MergeUpSpaceSAE(unsigned int k,
                                double space_best_score,
                                std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                                std::deque<std::shared_ptr<SAEFQIterator>> &fq_stack,
                                std::deque<unsigned int> &iterate_times,
                                SpaceHeap *collection)
{
  this->SplitIntoTwoSpace(k,space_best_score,record,fq_stack,iterate_times,collection);
  if(!fq_stack.empty()) {
    auto &back_fq = fq_stack.back();
    back_fq->MergeUpSpaceSAE(k, space_best_score, record, fq_stack, iterate_times, collection);
  }
}

void SAEFQIterator::GetFirst(unsigned int k) {
  double cost = 0;
  for (unsigned int j = 0; j < this->fr_ow_iterators_.size(); j++) {
    cost += this->fr_ow_iterators_[j]->top_element_ptr_->cost;
  }

  auto ele_ptr = std::unique_ptr<SAEElement>(new SAEElement);
  ele_ptr->cost = cost;
  this->top_element_ptr_ = std::move(ele_ptr);
}

void SAEFQIterator::Insert(std::shared_ptr<SAEList> FR_OW_iterator) {
  this->fr_ow_iterators_.push_back(std::move(FR_OW_iterator));
}

/**
 * Insert a bulk of FR or OW iterators.
 * inserting one certain type each time each time
 * certain type [i] specified by it's the i-th child of its father
 * @param FR_OW_iterators
 */
void SAEFQIterator::Insert(std::vector<std::shared_ptr<SAEList>> FR_OW_iterators) {
  this->fr_ow_iterators_ = std::move(FR_OW_iterators);
}

void SAEFQIterator::GetResult(std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
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
SAEFQIterator::SplitIntoTwoSpace(unsigned int k,
                                  double space_best_score,
                                  std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                                  std::deque<std::shared_ptr<SAEFQIterator>> &fq_stack,
                                  std::deque<unsigned int> &iterate_times,
                                  SpaceHeap *collection) {
  record.push_back(this->node_id_);
  fq_stack.push_back(this->shared_from_base<SAEFQIterator>());
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
    fr_ow_iterator->SplitIntoTwoSpace(k,space_best_score,record,fq_stack,iterate_times,collection);
    iterate_times.back()++;
  }
  iterate_times.pop_back();
  fq_stack.pop_back();
}

void SAEFQIterator::MergeUpSpaceSAE(unsigned int k,
                                    double space_best_score,
                                    std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                                    std::deque<std::shared_ptr<SAEFQIterator>> &fq_stack,
                                    std::deque<unsigned int> &iterate_times,
                                    SpaceHeap *collection) {

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
    this->fr_ow_iterators_[i]->SplitIntoTwoSpace(k,space_best_score,record,fq_stack,iterate_times,collection);
  }
  fq_stack.pop_back();
  iterate_times.pop_back();
  if(!fq_stack.empty())
  {
    auto upper_fq = fq_stack.back();
    upper_fq->MergeUpSpaceSAE(k,space_best_score,record,fq_stack,iterate_times,collection);
  }
}
