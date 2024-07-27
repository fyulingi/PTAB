
#ifndef PRIVATETOPK_QUERY_TOPK_DPP_GLOBALQUEUE_H_
#define PRIVATETOPK_QUERY_TOPK_DPP_GLOBALQUEUE_H_
#include "./DPPList.h"

class QueueElement{
 public:
  static const size_t InValidPos = -1;
  size_t pos_;
  enum class Type{FR,FQ,OW};
  // this element type, not parent type
  Type type_;
  /*
   * If it is FQ , then dpp_list_ is FQ itself
   * if is FR/OW , it is the FQ it belongs to
   */
  DPPList* fq_iterator_;
  DPPList* self_dpp_list_;
  // useful for FR and OW
  unsigned i_th_type;
  std::unique_ptr<element> i_element_;
  std::unique_ptr<FqElement> fq_element_;
  QueueElement() = default;
  double Score(){
    if(type_ == Type::FQ)
      return fq_element_->cost;
    else
      return i_element_->cost;
  }
};

// min top
class GlobalQueue {
 private:
  std::vector<std::unique_ptr<QueueElement>> contents_;
  void HeapSwap(int i, int j);
  void Down(size_t pos);
  void Up(size_t pos);
 public:
  void UpdateNode(std::unique_ptr<QueueElement> new_element,size_t pos);
  void Insert(std::unique_ptr<QueueElement> queue_element);
  size_t Size(){return contents_.size();}
  bool Empty(){return contents_.empty();}
  std::unique_ptr<QueueElement> ExtractMin();
  std::unique_ptr<QueueElement> Pop(size_t pos);
  GlobalQueue() = default;
};

#endif //PRIVATETOPK_QUERY_TOPK_DPP_GLOBALQUEUE_H_
