select ?s1 ?s2 ?s3 ?s4 ?o where{
    ?s1 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/OfficeHolder>.
    ?s2 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/University>.
    ?s3 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/Settlement>.
    ?s4 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/Country>.
    ?s1 <http://dbpedia.org/ontology/almaMater> ?s2.
    ?s2 <http://dbpedia.org/ontology/city> ?s3.
    ?s3 <http://dbpedia.org/ontology/country> ?s4.
    ?s4 <http://dbpedia.org/ontology/leaderName> ?s1.
    ?s2 <http://dbpedia.org/ontology/numberOfStudents> ?n.
    ?s4 <http://dbpedia.org/ontology/populationDensity> ?o.
} order by desc(?o+?n) limit 20