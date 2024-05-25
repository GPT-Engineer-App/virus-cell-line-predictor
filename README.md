# virus-cell-line-predictor

Build a tool that can predict a cell line suitable to culture a new virus using its genomic sequence data. Here is an example in Python: import pandas as pd
from Bio import Entrez, SeqIO
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score, classification_report
from sklearn.ensemble import RandomForestClassifier
from tqdm.notebook import tqdm
import pickle

# Set your email here
Entrez.email = "your_email@example.com"

# Define the list of viruses and their associated cell lines
virus_cell_lines = [
    {"virus": "Adenovirus", "cell_lines": ["A549", "BHK 21 (clone 13)", "HeLa", "MDCK", "MRC-5", "NCI-H292", "Wi 38"], "infectious_agents": ["Human adenovirus D", "Adenovirus type 3"]},
    {"virus": "HSV", "cell_lines": ["A549", "CV-1", "HeLa", "McCoy", "MRC-5", "NCI-H292", "Vero", "Vero76", "Wi 38"], "infectious_agents": ["HSV-1", "HSV-2"]},
    {"virus": "Influenza", "cell_lines": ["A549", "BHK 21 (clone 13)", "MDCK", "MRC-5", "NCI-H292", "Vero", "Wi 38"], "infectious_agents": ["Influenza A", "Influenza B"]},
    {"virus": "Measles", "cell_lines": ["A549", "CV-1", "NCI-H292", "Vero", "Vero/hSLAM"], "infectious_agents": ["Measles morbillivirus"]},
    {"virus": "Mumps", "cell_lines": ["A549", "MRC-5", "Vero", "Wi 38"], "infectious_agents": ["Mumps orthorubulavirus"]},
    {"virus": "Parainfluenza", "cell_lines": ["A549", "BHK 21 (clone 13)", "Vero"], "infectious_agents": ["Parainfluenza virus type 5"]},
    {"virus": "Poliovirus", "cell_lines": ["A549", "HeLa", "LLCMK2", "MRC-5", "Vero"], "infectious_agents": ["Poliovirus type I", "Poliovirus type 3"]},
    {"virus": "RSV", "cell_lines": ["A549", "CV-1", "MRC-5", "NCI-H292"], "infectious_agents": ["Human orthopneumovirus"]},
    {"virus": "Rotavirus", "cell_lines": ["A549", "Vero"], "infectious_agents": ["Rotavirus A"]},
    {"virus": "VZV", "cell_lines": ["A549", "CV-1", "HeLa", "MRC-5", "Wi 38"], "infectious_agents": ["Human alphaherpesvirus 3"]},
    {"virus": "MPV", "cell_lines": ["A549"], "infectious_agents": ["Human metapneumovirus"]},
    {"virus": "Reovirus", "cell_lines": ["BHK 21 (clone 13)", "MDCK", "NCI-H292"], "infectious_agents": ["Mammalian orthoreovirus"]},
    {"virus": "Rabies", "cell_lines": ["BHK 21 (clone 13)"], "infectious_agents": ["Rabies lyssavirus"]},
    {"virus": "Foot and Mouth", "cell_lines": ["BHK 21 (clone 13)"], "infectious_agents": ["Foot-and-mouth disease virus"]},
    {"virus": "Rubella", "cell_lines": ["BHK 21 (clone 13)", "Vero"], "infectious_agents": ["Rubella virus"]},
    {"virus": "CMV", "cell_lines": ["HeLa", "MRC-5", "Wi 38"], "infectious_agents": ["Human betaherpesvirus 5"]},
    {"virus": "Echovirus", "cell_lines": ["HeLa", "MRC-5", "Wi 38"], "infectious_agents": ["Echovirus"]},
    {"virus": "Rhinovirus", "cell_lines": ["HeLa", "LLCMK2", "MRC-5", "NCI-H292", "Wi 38"], "infectious_agents": ["Rhinovirus A"]},
    {"virus": "Enterovirus", "cell_lines": ["LLCMK2", "NCI-H292"], "infectious_agents": ["Enterovirus A"]},
    {"virus": "Poxvirus", "cell_lines": ["LLCMK2"], "infectious_agents": ["Orthopoxvirus"]},
    {"virus": "Vaccinia", "cell_lines": ["NCI-H292"], "infectious_agents": ["Vaccinia virus"]},
    {"virus": "BK polyomavirus", "cell_lines": ["NCI-H292"], "infectious_agents": ["BK polyomavirus"]},
    {"virus": "Coxsackie B", "cell_lines": ["Vero", "Vero76"], "infectious_agents": ["Coxsackievirus B"]},
    {"virus": "West Nile", "cell_lines": ["Vero76"], "infectious_agents": ["West Nile virus"]},
    {"virus": "Dengue", "cell_lines": ["BHK 21 (clone 13)"], "infectious_agents": ["Dengue virus"]},
    {"virus": "Vesicular stomatitis", "cell_lines": ["BHK 21 (clone 13)", "HeLa"], "infectious_agents": ["Vesicular stomatitis virus (Indiana strain)"]},
    {"virus": "Usutu", "cell_lines": ["Vero", "PK-15", "GEF"], "infectious_agents": ["Usutu virus"]},
    {"virus": "SARS-CoV-2", "cell_lines": ["Vero E6", "Caco-2/AT", "HuH-6/AT"], "infectious_agents": ["SARS-CoV-2"]},
    {"virus": "HIV", "cell_lines": ["CD4-positive lymphoid cell line", "HeLa"], "infectious_agents": ["Human immunodeficiency virus"]},
    {"virus": "HPV", "cell_lines": ["Various cell lines"], "infectious_agents": ["Human papillomavirus type 18"]},
    {"virus": "XMLV", "cell_lines": ["Various cell lines"], "infectious_agents": ["Xenotropic murine leukemia virus"]},
    {"virus": "Bovine polyoma", "cell_lines": ["SK-BR-3"], "infectious_agents": ["Bovine polyomavirus"]},
    {"virus": "Feline sarcoma", "cell_lines": ["Various cell lines"], "infectious_agents": ["Feline sarcoma virus"]},
    {"virus": "Mason-Pfizer monkey", "cell_lines": ["Various cell lines"], "infectious_agents": ["Mason-Pfizer monkey virus"]},
    {"virus": "Squirrel monkey retrovirus", "cell_lines": ["Various cell lines"], "infectious_agents": ["Squirrel monkey retrovirus"]}
]

# Create a DataFrame
data = {
    "Cell Line": [],
    "Infectious Agents": [],
    "Virus": []
}

for entry in virus_cell_lines:
    virus = entry["virus"]
    for cell_line in entry["cell_lines"]:
        for agent in entry["infectious_agents"]:
            data["Cell Line"].append(cell_line)
            data["Infectious Agents"].append(agent)
            data["Virus"].append(virus)

df = pd.DataFrame(data)

# Function to fetch genomic sequence data from NCBI with retries
def fetch_genomic_sequence(virus_name, retries=3):
    for attempt in range(retries):
        try:
            handle = Entrez.esearch(db="nucleotide", term=f"{virus_name}[Organism] AND complete genome", retmax=1)
            record = Entrez.read(handle)
            handle.close()
            if record["IdList"]:
                seq_id = record["IdList"][0]
                handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
                seq_record = SeqIO.read(handle, "fasta")
                handle.close()
                return str(seq_record.seq)
            else:
                return None
        except Exception as e:
            print(f"Error fetching data for {virus_name} (attempt {attempt + 1}): {e}")
            time.sleep(5)  # Wait before retrying
    return None

# Fetch genomic sequences for all viruses
print("Fetching Genomic Sequences...")
df["Genomic Sequence"] = [fetch_genomic_sequence(agent) for agent in tqdm(df["Infectious Agents"], desc="Fetching Genomic Sequences")]

# Check how many sequences were fetched
num_sequences_fetched = df["Genomic Sequence"].notna().sum()
print(f"Number of sequences fetched: {num_sequences_fetched}")

# Drop rows with missing genomic sequences
df.dropna(subset=["Genomic Sequence"], inplace=True)

# Feature extraction: For simplicity, we use GC content as a feature
df['GC Content'] = df['Genomic Sequence'].apply(lambda seq: (seq.count('G') + seq.count('C')) / len(seq))

# Encode the cell line names
label_encoder = LabelEncoder()
df["Cell Line Encoded"] = label_encoder.fit_transform(df["Cell Line"])

# Remove classes with fewer than 2 instances
class_counts = df["Cell Line Encoded"].value_counts()
df = df[df["Cell Line Encoded"].isin(class_counts[class_counts >= 2].index)]

# Prepare features and labels
X = df[['GC Content']]
y = df['Cell Line Encoded']

# Split the data into training and testing sets
print("Splitting Data...")
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

# Initialize and train the Random Forest classifier
clf = RandomForestClassifier(n_estimators=100, random_state=42)
clf.fit(X_train, y_train)

# Make predictions on the test set
y_pred = clf.predict(X_test)

# Evaluate the model
accuracy = accuracy_score(y_test, y_pred)
print(f'Accuracy: {accuracy * 100:.2f}%')

# Classification report
print("Classification Report:")
print(classification_report(y_test, y_pred, target_names=label_encoder.classes_, labels=np.unique(y_test)))

# Function to predict cell line for a new sequence using GC content
def predict_cell_line(sequence):
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
    prediction = clf.predict([[gc_content]])
    return label_encoder.inverse_transform(prediction)[0]

# Example usage
new_sequence = "ATGCGTACGTAGCTAGCTAGCTA"
predicted_cell_line = predict_cell_line(new_sequence)
print(f'Predicted cell line for the new sequence: {predicted_cell_line}')

# Save the model
with open('cell_line_prediction_model.pkl', 'wb') as handle:
    pickle.dump(clf, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Save the label encoder
with open('label_encoder.pickle', 'wb') as handle:
    pickle.dump(label_encoder, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Save the DataFrame to a CSV file
df.to_csv("virus_cell_line_data.csv", index=False)

# Display the DataFrame
print(df.head())

# t-SNE visualization
print("Performing t-SNE visualization...")
perplexity = min(30, len(X_train) - 1)  # Ensure perplexity is less than the number of samples
X_embedded_tsne = TSNE(n_components=2, n_iter=2000, init='random', perplexity=perplexity).fit_transform(X_train)
plt.figure(figsize=(10, 6))
scatter = plt.scatter(X_embedded_tsne[:, 0], X_embedded_tsne[:, 1], c=np.argmax(y_train, axis=1), cmap='viridis', s=5)
plt.colorbar(scatter, label='Cell Line')
plt.title('t-SNE visualization of genomic sequences')
plt.xlabel('t-SNE component 1')
plt.ylabel('t-SNE component 2')
plt.show()

# Save the t-SNE plot
plt.savefig("tsne_visualization.png")

# Print model summary
print("Model Summary:")
model.summary()

# Classification report
y_pred = model.predict(X_test)
y_pred_classes = np.argmax(y_pred, axis=1)
y_test_classes = np.argmax(y_test, axis=1)
print("Classification Report:")
print(classification_report(y_test_classes, y_pred_classes, target_names=label_encoder.classes_, labels=np.unique(y_test_classes)))


## Collaborate with GPT Engineer

This is a [gptengineer.app](https://gptengineer.app)-synced repository 🌟🤖

Changes made via gptengineer.app will be committed to this repo.

If you clone this repo and push changes, you will have them reflected in the GPT Engineer UI.

## Tech stack

This project is built with React and Chakra UI.

- Vite
- React
- Chakra UI

## Setup

```sh
git clone https://github.com/GPT-Engineer-App/virus-cell-line-predictor.git
cd virus-cell-line-predictor
npm i
```

```sh
npm run dev
```

This will run a dev server with auto reloading and an instant preview.

## Requirements

- Node.js & npm - [install with nvm](https://github.com/nvm-sh/nvm#installing-and-updating)
