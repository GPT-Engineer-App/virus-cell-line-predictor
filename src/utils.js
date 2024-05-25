export const fetchGenomicSequence = async (virusName, retries = 3) => {
  const apiUrl = `https://api.ncbi.nlm.nih.gov/genomes?term=${virusName}[Organism]&retmax=1`;
  for (let attempt = 0; attempt < retries; attempt++) {
    try {
      const response = await fetch(apiUrl);
      const data = await response.json();
      if (data && data.length > 0) {
        return data[0].sequence;
      } else {
        return null;
      }
    } catch (e) {
      console.error(`Error fetching data for ${virusName} (attempt ${attempt + 1}): ${e}`);
      await new Promise((resolve) => setTimeout(resolve, 5000));
    }
  }
  return null;
};

export const predictCellLine = async (sequence) => {
  const gcContent = (sequence.match(/G/g).length + sequence.match(/C/g).length) / sequence.length;

  if (gcContent > 0.5) {
    return "High GC Content Cell Line";
  } else {
    return "Low GC Content Cell Line";
  }
};
