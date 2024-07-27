select ?s1 ?s2 ?s3 ?s4 ?o where{
    ?s1 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/Currency>.
    ?s2 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/Country>.
    ?s3 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/OfficeHolder>.
    ?s4 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/Country>.
    ?s1 <http://dbpedia.org/ontology/usingCountry> ?s2.
    ?s2 <http://dbpedia.org/ontology/leaderName> ?s3.
    ?s3 <http://dbpedia.org/ontology/birthPlace> ?s4.
    ?s4 <http://dbpedia.org/ontology/currency> ?s1.
    ?s2 <http://dbpedia.org/ontology/areaTotal> ?n.
    ?s3 <http://dbpedia.org/property/percentage> ?m.
    ?s4 <http://dbpedia.org/ontology/areaTotal> ?o.
} order by desc(?n+?m+?o) limit 20