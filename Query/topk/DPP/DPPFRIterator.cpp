
#include "DPPFQIterator.h"
#include "DPPFRIterator.h"
#include "DPPOWIterator.h"


// getFirst() also implemented here
void DPPFRIterator::TryGetNext(unsigned int k,GlobalQueue *global_queue,
                               DPPFQIterator* fq_iterator,unsigned int i_th_child) {
  // get first is not necessary here
  auto m = this->pool_.size();
  auto em = this->pool_.back();
  auto em_pointer = this->fqs_map_[em.node];
  this->queue_.popMin();
  if(NextEPoolElement(k,em_pointer,em.index+1,global_queue))
  {
    decltype(em) e;
    e.node = em.node;
    e.index = em.index + 1;
    e.cost = em.cost + DeltaCost(em_pointer,e.index);
    this->queue_.push(e);
  }
  if(!queue_.empty()) {
    if(this->queue_.size()>k-m)
      this->queue_.popMax();
  }
  if(this->queue_.empty()) {
    if(this->GQ_pos_!=QueueElement::InValidPos)
      global_queue->Pop(this->GQ_pos_);
  }
  else
  {
    auto e = this->queue_.findMin();
    auto queue_element = std::unique_ptr<QueueElement>(new QueueElement);
    queue_element->fq_iterator_ = (DPPList*)fq_iterator;
    queue_element->self_dpp_list_ = this;
    queue_element->i_element_ = std::unique_ptr<element>(new element(e));
    queue_element->i_th_type = i_th_child;
    queue_element->type_ = QueueElement::Type::FR;
    global_queue->Insert(std::move(queue_element));
  }
}


/**
 * Insert One FQ iterator into the FR
 * @param k
 * @param fq_pointer
 */
void DPPFRIterator::Insert(TYPE_ENTITY_LITERAL_ID fq_id,
                           std::shared_ptr<DPPFQIterator> fq_pointer,
                           OnePointPredicatePtr predicates_vec) {
  (*this->type_predicates_)[fq_id] = predicates_vec;
  auto cost = fq_pointer->pool_[0].cost;
  element e{};
  e.cost = cost;
  e.index = 0;
  e.node = fq_id;
  this->fqs_map_[fq_id] = fq_pointer;
  queue_.push(e);
}

void DPPFRIterator::Register(TYPE_ENTITY_LITERAL_ID fq_id,
                           std::shared_ptr<DPPFQIterator> fq_pointer,
                           OnePointPredicatePtr predicates_vec) {
  (*this->type_predicates_)[fq_id] = predicates_vec;
  this->fqs_map_[fq_id] = fq_pointer;
}

double DPPFRIterator::DeltaCost(std::shared_ptr<DPPFQIterator> node_pointer, int index) {
  auto delta =  node_pointer->pool_[index].cost - node_pointer->pool_[index-1].cost;
  return delta;
}

bool DPPFRIterator::NextEPoolElement(unsigned int k,
                                     std::shared_ptr<DPPFQIterator> node_pointer,
                                     unsigned int index,
                                     GlobalQueue *global_queue) {
  if(index< node_pointer->pool_.size())
    return true ;

  if(index == node_pointer->pool_.size()) {
    // line 14 in algo.7
    if(node_pointer->TempParents().empty())
      node_pointer->TryGetNext(k,global_queue);
    // line 15 in algo.7
    // but different from the paper, the index should be (iPoolJ size - 1)
    // as in line 6 in algo.7 it reads iPoolJ[index], so the index should
    // represent the last element in the pool
    node_pointer->TempParents().emplace_back(this,this->pool_.size()-1);
  }
  return false;
}

/**
 *
 * @param i_th
 * @param record
 * @param predicate_information not used, because FR itself saves the information
 */
void DPPFRIterator::GetResult(int i_th, std::shared_ptr<std::vector<TYPE_ENTITY_LITERAL_ID>> record,
                              NodeOneChildVarPredicatesPtr predicate_information) {
  auto fq = this->pool_[i_th];
  auto fq_i_th = fq.index;
  auto fq_id = fq.node;
  auto fq_pointer = this->fqs_map_[fq_id];
  if(this->type_predicates_->find(fq_id)!=this->type_predicates_->end()) {
    auto predicates = (*this->type_predicates_)[fq_id];
    for (auto pre_id:*predicates)
      record->push_back(pre_id);
  }
  fq_pointer->GetResult(fq_i_th,record);
}


DPPFRIterator::DPPFRIterator() {
  this->type_predicates_=std::make_shared<NodeOneChildVarPredicates>();
}

void DPPFRIterator::dInsert(unsigned int k,
                            DPPFQIterator* child_fq_pointer,
                            unsigned i_th_type,
                            GlobalQueue *global_queue) {
  // line 1-4 in algo 7.
  auto cost = child_fq_pointer->pool_[0].cost;
  TYPE_ENTITY_LITERAL_ID fq_id = child_fq_pointer->GetID();
  element e{};
  e.cost = cost;
  e.index = 0;
  e.node = fq_id;
  queue_.push(e);

  if(this->queue_.size()> k)
    this->queue_.popMax();

  GQ_Update(global_queue, e, i_th_type);
}

void DPPFRIterator::rInsert(unsigned int k,
                            int index,
                            DPPList* fq_pointer,
                            unsigned i_th_type,
                            GlobalQueue *global_queue) {
  // line 5-7 in algo.7
  auto e_old = this->pool_[index];
  element e{};
  e.index = e_old.index+1;
  e.cost = e_old.cost + DeltaCost(this->fqs_map_[ e_old.node],e.index);
  e.node = e_old.node;
  queue_.push(e);

  if(this->queue_.size()> k)
    this->queue_.popMax();

  GQ_Update(global_queue,e,i_th_type);
}

void DPPFRIterator::GQ_Update(GlobalQueue *global_queue,element e, unsigned i_th_type) {
  if(e != this->queue_.findMin()) return;
  auto queue_element = std::unique_ptr<QueueElement>(new QueueElement);
  queue_element->fq_iterator_ = this->parent_fq_;
  queue_element->self_dpp_list_ = this;
  queue_element->i_element_ = std::unique_ptr<element>(new element);
  queue_element->i_element_->cost = e.cost;
  queue_element->i_element_->node = e.node;
  queue_element->i_element_->index = e.index;
  queue_element->i_th_type = i_th_type;
  queue_element->type_ = QueueElement::Type::FR;
  // line 25 in algo.5
  if(this->GQ_pos_ == QueueElement::InValidPos)
  {
    global_queue->Insert(std::move(queue_element));
  }
  else
  {
    global_queue->UpdateNode(std::move(queue_element),this->GQ_pos_);
  }
  // line 26 in algo.5 is in the global_queue->Insert
}
void DPPFRIterator::GetFirst() {
  if(this->queue_.empty())
    return;
  if(this->pool_.size()==0)
  {
    this->pool_.push_back(this->queue_.findMin());
    return;
  }
}

