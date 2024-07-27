select ?s1 ?s2 ?s3 ?s4 ?o where{
    ?s1 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/SoccerClub>.
    ?s2 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/SoccerManager>.
    ?s3 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/SoccerClub>.
    ?s4 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/Stadium>.
    ?s1 <http://dbpedia.org/ontology/manager> ?s2.
    ?s2 <http://dbpedia.org/ontology/managerClub> ?s3.
    ?s3 <http://dbpedia.org/ontology/ground> ?s4.
    ?s4 <http://dbpedia.org/ontology/tenant> ?s1.
    ?s1 <http://dbpedia.org/ontology/capacity> ?n.
    ?s3 <http://dbpedia.org/ontology/capacity> ?m.
    ?s4 <http://www.w3.org/2003/01/geo/wgs84_pos#lat> ?o.
} order by desc(?n-?m+?o) limit 20