select distinct ?e ?t1
where{
   ?e <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> ?t1.
   ?e <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> ?t2.
   filter(?t1 != ?t2)
} order by ?e