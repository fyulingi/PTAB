
#include "SpaceSAE.h"

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
  auto fq_pointer = this->top_element_ptr_->fq_;
  if(this->top_element_ptr_->predicate_!= nullptr){
    for (auto pre_id:*this->top_element_ptr_->predicate_)
      record->push_back(pre_id);
  }
  fq_pointer->GetResult(record);
}

/**
 * if FR = {a1 a2 a3 ... an}
 * then split it into
 * { a1 } and { a2 a3 ... an }
 */
void
SAEFRIterator::SplitIntoTwoSpace(unsigned int k,
                                 double space_best_score,
                                 std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                                 std::unique_ptr<std::stack<std::shared_ptr<SAEList>>> fr_ow_iterators,
                                 unsigned int iterate_time,
                                 SpaceHeapSAE *collection) {
  auto top_id = this->top_element_ptr_->node;
  auto top_cost_ = this->top_element_ptr_->cost;
  // other node + best score of  { a2 a3 ... an } =  remaining_cost + best of  { a2 a3 ... an }
  auto remaining_cost = space_best_score - top_cost_;

  std::shared_ptr<SAEFRIterator> the_other_fr;
  if(this->next_!=nullptr){
    the_other_fr = this->next_;
  }
  else {
    the_other_fr = std::make_shared<SAEFRIterator>();
    for (unsigned int i = 0; i < this->fq_ids_.size(); i++){
      if(this->fq_ids_[i] != top_id) {
        the_other_fr->Insert(k, this->fq_ids_[i],
                             this->fqs_[i],
                             this->predicates_[i]);
      }
    }
  }
  if(this->fq_ids_.size()>1) // the other_fr != empty
  {
    fr_ow_iterators->pop();
    fr_ow_iterators->push(the_other_fr);
    auto new_sub_space = std::unique_ptr<SpaceSAE>(new SpaceSAE(record,
                                                                iterate_time,
                                                                std::move(fr_ow_iterators),
                                                                remaining_cost + the_other_fr->top_element_ptr_->cost));
    collection->push(std::move(new_sub_space));
  }
  if(this->next_==nullptr)
    this->next_ = the_other_fr;
  this->SplitTopFq(k,space_best_score,
                     record,collection);
}

SAEFRIterator::SAEFRIterator() {
}

void SAEFRIterator::SplitTopFq(unsigned int k,
                               double space_best_score,
                               std::vector<TYPE_ENTITY_LITERAL_ID> &record,
                               SpaceHeapSAE *collection){
  auto fq = this->top_element_ptr_->fq_;
  if(this->top_element_ptr_->predicate_!= nullptr){
    for (auto pre_id:*this->top_element_ptr_->predicate_)
      record.push_back(pre_id);
  }
  std::unique_ptr<std::stack<std::shared_ptr<SAEList>>> the_only_fq(
      new std::stack<std::shared_ptr<SAEList>>()
      );
  the_only_fq->push(fq);
  fq->SplitIntoTwoSpace(k,space_best_score,record,std::move(the_only_fq),0,collection);
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
  SAEElement e{};
  e.cost = score;
  e.node = id;
  this->pool_.push_back(e);
  if(this->pool_.size()==1)
    this->top_element_ptr_ = std::unique_ptr<SAEElement>(new SAEElement(this->pool_[0]));
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
                                 std::vector<TYPE_ENTITY_LITERAL_ID>& record,
                                 std::unique_ptr<std::stack<std::shared_ptr<SAEList>>> fr_ow_iterators,
                                 unsigned int iterate_time,
                                 SpaceHeapSAE *collection) {
  auto top_cost_ = this->pool_[0].cost;
  auto remaining_cost = space_best_score - top_cost_;
  auto the_other_ow = std::make_shared<SAEOWIterator>();
  bool the_other_empty = this->pool_.size() == 1;

  if(this->next_!=nullptr)
    the_other_ow = this->next_;
  else
  {
    for (unsigned int i = 1; i < this->pool_.size(); i++)
      the_other_ow->RankedInsert(k, this->pool_[i].node, this->pool_[i].cost);
  }
  if(!the_other_empty)
  {
    auto new_sub_space = std::unique_ptr<SpaceSAE>(new SpaceSAE(record,
                                                                iterate_time,
                                                                std::move(fr_ow_iterators),
                                                                remaining_cost + the_other_ow->top_element_ptr_->cost));
    collection->push(std::move(new_sub_space));
  }
  if(this->next_==nullptr)
    this->next_ = the_other_ow;

  record.push_back(this->pool_[0].node);

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


/**
 * 对一个FQ的空间进行分割，要完成两部分事情:
 * 1. 收集结果
 * 2. 把深度优先搜索的过程传下去
 * 3. 生成一个新的FQ（和一个包裹的FR），插入collection中
 * @param k
 * @param space_best_score
 * @param record
 * @param fr_ow_iterators contains the shared pointer of this fq
 * @param iterate_time ignored in this function
 * @param collection
 */
void
SAEFQIterator::SplitIntoTwoSpace(unsigned int k,
                                 double space_best_score,
                                 std::vector<TYPE_ENTITY_LITERAL_ID>& record_up_layer,
                                 std::unique_ptr<std::stack<std::shared_ptr<SAEList>>> fr_ow_iterators,
                                 unsigned int iterate_time,
                                 SpaceHeapSAE *collection) {
  auto remaining_cost = space_best_score - this->node_score_;
  record_up_layer.push_back(this->node_id_);
  iterate_time = 0;
  std::vector<TYPE_ENTITY_LITERAL_ID> record;
  // 这个地方不能调用  this->fr_ow_iterators_ 的split，而应该调用
  // 传入的 std::unique_ptr<std::stack<std::shared_ptr<SAEList>>>  fr_ow_iterators
  auto InsertSpace = [&record, k, this, space_best_score]
      (std::unique_ptr<std::stack<std::shared_ptr<SAEList>>> &fr_ow_iterators,
       int pos,
       SpaceHeapSAE *collection) {
    auto fr_ow_iterator = fr_ow_iterators->top();
    if (fr_ow_iterator->Type() == OrderedListType::OW) {
      auto ow_predicates = this->types_predicates_[pos];
      auto node_id = fr_ow_iterator->top_element_ptr_->node;
      auto predicates = (*ow_predicates)[node_id];
      for (auto predicate_id:*predicates)
        record.push_back(predicate_id);
    }
    std::unique_ptr<std::stack<std::shared_ptr<SAEList>>>
        unique(new std::stack<std::shared_ptr<SAEList>>(*fr_ow_iterators));
    fr_ow_iterator->SplitIntoTwoSpace(k, space_best_score, record, std::move(unique), pos, this->small_queue_);
  };

  // 如果之前已经分裂过了，拿best solution继续分裂即可,
  // 这个新空间插入的是这个FQ的优先队列
  if(this->small_queue_ != nullptr && !this->small_queue_->empty())
  {
    auto best_solution = this->small_queue_->popMin();
    for(auto solution_point: best_solution->one_node_space_)
      record_up_layer.push_back(solution_point);
    unsigned int iterated_size = this->fr_ow_iterators_.size() -  best_solution->fr_ow_iterators_->size();
    for(unsigned int i = 0; !best_solution->fr_ow_iterators_->empty();i++)
    {
      InsertSpace(best_solution->fr_ow_iterators_,
                  i + iterated_size,
                  this->small_queue_);
      best_solution->fr_ow_iterators_->pop();
    }
  } // 若之前没分配过，则需要初始化，并且插入
  else if(this->small_queue_ == nullptr){
    this->small_queue_ = new SpaceHeapSAE;
    // 初始化 small_queue_
    std::unique_ptr<std::stack<std::shared_ptr<SAEList>>>
        unique(new std::stack<std::shared_ptr<SAEList>>());
    for(int j = this->fr_ow_iterators_.size()-1;j>=0;j--)
      unique->push(this->fr_ow_iterators_[j]);
    for(unsigned int i=0;i<this->fr_ow_iterators_.size();i++)
    {
      InsertSpace(unique,i,
                  this->small_queue_);
      unique->pop();
    }
  }
  for(auto solution_point: record)
    record_up_layer.push_back(solution_point);
  if(this->small_queue_ != nullptr && !this->small_queue_->empty()) {
    auto &min_ele = this->small_queue_->findMin();
    auto new_sub_space = std::unique_ptr<SpaceSAE>(new SpaceSAE(record_up_layer,
                                                                0,
                                                                std::move(fr_ow_iterators),
                                                                min_ele->best_score_
    ));
    collection->push(std::move(new_sub_space));
  }
}


void SpaceSAE::Set(std::vector<TYPE_ENTITY_LITERAL_ID> one_node_space,
                   unsigned int iterated_items,
                   std::unique_ptr<std::stack<std::shared_ptr<SAEList>>> fr_ow_iterators,
                   double best_score) {
  this->one_node_space_ = std::move(one_node_space);
  this->iterated_items_ = iterated_items;
  this->fr_ow_iterators_ = std::move(fr_ow_iterators);
  this->best_score_ = best_score;
}
