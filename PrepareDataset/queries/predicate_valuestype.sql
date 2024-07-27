select distinct ?p ?time where
{
    select ?p (count(?p) as ?time) where{
        ?s ?p ?o.
        filter( datatype(?o)=<http://www.w3.org/2001/XMLSchema#double> || datatype(?o)=<http://www.w3.org/2001/XMLSchema#integer> || datatype(?o)=<http://www.w3.org/2001/XMLSchema#float> || datatype(?o)=<http://dbpedia.org/datatype/usDollar> || datatype(?o)=<http://www.w3.org/2001/XMLSchema#positiveInteger> || datatype(?o)=<http://www.w3.org/2001/XMLSchema#nonNegativeInteger> || datatype(?o)=<http://dbpedia.org/datatype/euro> ) 
    } group by ?p
}order by desc(?time)
