select *
where
{
        select  ?o  (count(*) as ?time) 
        where
        {
            ?s  <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> ?o.
        }group by ?o
}order by desc(?time)