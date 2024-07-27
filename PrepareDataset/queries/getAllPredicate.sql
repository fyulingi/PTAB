select *
where
{
        select  ?p  (count(*) as ?time) 
        where
        {
            ?s ?p ?o.
        }group by ?p
}order by desc(?time)