PREFIX  dc:   <http://purl.org/dc/elements/1.1/>
PREFIX  dbpedia2: <http://dbpedia.org/property/>
PREFIX  rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX  foaf: <http://xmlns.com/foaf/0.1/>
PREFIX  dbpedia: <http://dbpedia.org/>
PREFIX  owl:  <http://www.w3.org/2002/07/owl#>
PREFIX  xsd:  <http://www.w3.org/2001/XMLSchema#>
PREFIX  rdf:  <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX  do:   <http://dbpedia.org/ontology/>
PREFIX  skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX  dr:   <http://dbpedia.org/resource/>

SELECT DISTINCT  ?a ?l ?w
WHERE
  { ?a rdf:type foaf:Person .
    ?a rdfs:label ?l.
    ?a <http://dbpedia.org/property/weight> ?w.
  } order by desc(?w)
LIMIT   5