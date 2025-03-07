PREFIX  dc:   <http://purl.org/dc/elements/1.1/>
PREFIX  :     <http://dbpedia.org/resource/>
PREFIX  rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX  dbpedia2: <http://dbpedia.org/property/>
PREFIX  foaf: <http://xmlns.com/foaf/0.1/>
PREFIX  owl:  <http://www.w3.org/2002/07/owl#>
PREFIX  xsd:  <http://www.w3.org/2001/XMLSchema#>
PREFIX  dbpedia: <http://dbpedia.org/>
PREFIX  rdf:  <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX  skos: <http://www.w3.org/2004/02/skos/core#>

SELECT DISTINCT  ?planet ?per ?apo 
WHERE
  { 
      ?planet rdf:type <http://dbpedia.org/ontology/Planet>.
      ?planet <http://dbpedia.org/ontology/periapsis> ?per.
       ?planet <http://dbpedia.org/ontology/apoapsis> ?apo.
    
  } order by desc(?per + ?apo)
limit 5