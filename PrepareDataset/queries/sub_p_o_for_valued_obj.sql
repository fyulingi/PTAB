select ?stype ?p ?time
where {
    select ?stype ?p ( count(*) as ?time) where
    {
        ?s <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> ?stype.
        filter( ?stype != <http://www.w3.org/2002/07/owl#Thing> )
    ?s ?p ?o.
    filter( datatype(?o)=<http://www.w3.org/2001/XMLSchema#double> || datatype(?o)=<http://www.w3.org/2001/XMLSchema#integer> || datatype(?o)=<http://www.w3.org/2001/XMLSchema#float> || datatype(?o)=<http://dbpedia.org/datatype/usDollar> || datatype(?o)=<http://www.w3.org/2001/XMLSchema#positiveInteger> || datatype(?o)=<http://www.w3.org/2001/XMLSchema#nonNegativeInteger> || datatype(?o)=<http://dbpedia.org/datatype/euro> )
    }group by ?stype ?p 
}order by desc(?time)