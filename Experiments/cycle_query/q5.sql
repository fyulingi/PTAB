select ?s1 ?s2 ?s3 ?s4 ?s5 ?o where{
    ?s1 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/City>.
    ?s2 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/Settlement>.
    ?s3 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/City>.
    ?s4 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/Settlement>.
    ?s5 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/Country>.
    ?s1 <http://dbpedia.org/ontology/part> ?s2.
    ?s2 <http://dbpedia.org/ontology/isPartOf> ?s3.
    ?s3 <http://dbpedia.org/ontology/part> ?s4.
    ?s4 <http://dbpedia.org/ontology/country> ?s5.
    ?s5 <http://dbpedia.org/ontology/capital> ?s1.
    ?s1 <http://dbpedia.org/property/areaKm> ?n.
    ?s5 <http://www.w3.org/2003/01/geo/wgs84_pos#lat> ?o.
} order by desc(?n-?o) limit 20
