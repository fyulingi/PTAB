PREFIX  dc:   <http://purl.org/dc/elements/1.1/>
PREFIX  :     <http://dbpedia.org/resource/>
PREFIX  rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX  dbpedia2: <http://dbpedia.org/property/>
PREFIX  foaf: <http://xmlns.com/foaf/0.1/>
PREFIX  dbo:  <http://dbpedia.org/ontology/>
PREFIX  owl:  <http://www.w3.org/2002/07/owl#>
PREFIX  xsd:  <http://www.w3.org/2001/XMLSchema#>
PREFIX  dbpedia: <http://dbpedia.org/>
PREFIX  rdf:  <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX  skos: <http://www.w3.org/2004/02/skos/core#>

SELECT  *
WHERE
  { 
    ?artist rdf:type dbo:MusicalArtist .
    ?artist dbo:birthPlace ?place .
    ?place rdf:type dbo:City .
    ?place <http://dbpedia.org/ontology/areaTotal> ?area.
  }order by Desc(?area)
LIMIT 5