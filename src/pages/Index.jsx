import React, { useState } from "react";
import { Container, VStack, Text, Input, Button, Box, Alert, AlertIcon } from "@chakra-ui/react";
import { FaVirus, FaSearch } from "react-icons/fa";
import { IconButton } from "@chakra-ui/react";

const Index = () => {
  const [virusName, setVirusName] = useState("");
  const [sequence, setSequence] = useState("");
  const [loading, setLoading] = useState(false);
  const [predictedCellLine, setPredictedCellLine] = useState("");
  const [error, setError] = useState("");

  const handleFetchSequence = async () => {
    setLoading(true);
    setError("");
    try {
      const seq = await fetchGenomicSequence(virusName);
      if (seq) {
        setSequence(seq);
      } else {
        setError("No sequence found for the given virus.");
      }
    } catch (err) {
      setError("Error fetching sequence data.");
    }
    setLoading(false);
  };

  const handlePredictCellLine = async () => {
    setLoading(true);
    setError("");
    try {
      const cellLine = await predictCellLine(sequence);
      setPredictedCellLine(cellLine);
    } catch (err) {
      setError("Error predicting cell line.");
    }
    setLoading(false);
  };

  return (
    <Container centerContent maxW="container.md" height="100vh" display="flex" flexDirection="column" justifyContent="center" alignItems="center">
      <VStack spacing={4}>
        <Text fontSize="2xl">Virus Cell Line Predictor</Text>
        <Box width="100%">
          <Input placeholder="Enter virus name" value={virusName} onChange={(e) => setVirusName(e.target.value)} />
          <IconButton aria-label="Fetch Sequence" icon={<FaSearch />} onClick={handleFetchSequence} isLoading={loading} mt={2} />
        </Box>
        {sequence && (
          <Box width="100%">
            <Text>Genomic Sequence:</Text>
            <Text fontSize="sm" wordBreak="break-all">
              {sequence}
            </Text>
            <Button leftIcon={<FaVirus />} onClick={handlePredictCellLine} isLoading={loading} mt={2}>
              Predict Cell Line
            </Button>
          </Box>
        )}
        {predictedCellLine && (
          <Alert status="success" mt={4}>
            <AlertIcon />
            Predicted Cell Line: {predictedCellLine}
          </Alert>
        )}
        {error && (
          <Alert status="error" mt={4}>
            <AlertIcon />
            {error}
          </Alert>
        )}
      </VStack>
    </Container>
  );
};

export default Index;

// utils.js (Assuming this file exists and contains the following functions)

// Set your email here
Entrez.email = "your_email@example.com";

// Function to fetch genomic sequence data from NCBI with retries
export const fetchGenomicSequence = async (virusName, retries = 3) => {
  for (let attempt = 0; attempt < retries; attempt++) {
    try {
      const handle = await Entrez.esearch({ db: "nucleotide", term: `${virusName}[Organism] AND complete genome`, retmax: 1 });
      const record = await Entrez.read(handle);
      handle.close();
      if (record.IdList.length > 0) {
        const seqId = record.IdList[0];
        const handle = await Entrez.efetch({ db: "nucleotide", id: seqId, rettype: "fasta", retmode: "text" });
        const seqRecord = await SeqIO.read(handle, "fasta");
        handle.close();
        return seqRecord.seq.toString();
      } else {
        return null;
      }
    } catch (e) {
      console.error(`Error fetching data for ${virusName} (attempt ${attempt + 1}): ${e}`);
      await new Promise((resolve) => setTimeout(resolve, 5000)); // Wait before retrying
    }
  }
  return null;
};

// Function to predict cell line for a new sequence using GC content
export const predictCellLine = async (sequence) => {
  const gcContent = (sequence.match(/G/g).length + sequence.match(/C/g).length) / sequence.length;
  const model = await pickle.load(open("cell_line_prediction_model.pkl", "rb"));
  const labelEncoder = await pickle.load(open("label_encoder.pickle", "rb"));
  const prediction = model.predict([[gcContent]]);
  return labelEncoder.inverse_transform(prediction)[0];
};
