select ?s1 ?s2 ?s3 ?s4 ?o where{
    ?s1 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/OfficeHolder>.
    ?s2 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/OfficeHolder>.
    ?s3 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/PoliticalParty>.
    ?s4 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/Country>.
    ?s1 <http://dbpedia.org/ontology/president> ?s2.
    ?s2 <http://dbpedia.org/ontology/party> ?s3.
    ?s3 <http://dbpedia.org/ontology/country> ?s4.
    ?s4 <http://dbpedia.org/ontology/leaderName> ?s1.
    ?s2 <http://www.w3.org/2003/01/geo/wgs84_pos#lat> ?n.
    ?s4 <http://dbpedia.org/ontology/areaTotal> ?o.
} order by desc(2*?n+?o) limit 20