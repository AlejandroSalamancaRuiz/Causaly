
# Gene-Disease Analysis AI Assistant

This AI Assistant is designed to help users understand how genes are involved in diseases. It has three main functions:

### Key Features:
1. **Retrieve Diseases with Specific Genes**: Identify diseases associated with particular genes from a curated database.
2. **Generate Hypotheses on Gene Roles in Diseases**: Analyze how individual genes interact within biological pathways, predict downstream gene effects, and hypothesize their potential contributions to diseases.
3. **Process Queries Involving Multiple Genes**: Provide insights into how a group of genes may collectively impact biological processes and contribute to disease outcomes.

### System Overview:
This AI Assistant is built using the LangChain framework, with a series of tools corresponding to each key function, making it a versatile platform for gene-disease analysis.

### How to Use:

To get started, follow these steps:

#### 1. Create a Data Folder
- In the root directory of the project, create a folder named `data`.

#### 2. Prepare KEGG Data
- Inside the `data` folder, create a folder named `KEGG_data`.
- Download the KEGG disease data and place the following subfolders inside `KEGG_data`:
  - `KGML`: This should contain the XML files with disease information.
  - `PNG`: This should contain pathway images (optional, if needed for visualization).

#### 3. Prepare Gene Ontology Data
- Inside the `data` folder, create another folder named `GOA_human`.
- Download the **Gene Association File** for *Homo sapiens* and place it in `GOA_human`.
- Download the **GO Basic Ontology** file (`go-basic.obo`) and place it in the `GOA_human` folder.

#### 4. Configure Environment Variables
- In the root directory, locate the `.env` file or create one if not already present.
- Set your OpenAI API key by modifying the `OPENAI_API_KEY` variable:

```bash
OPENAI_API_KEY=<your-openai-key>
```

#### 5. Process Data
- Run the `Data_Processing.py` file to process the KEGG and Gene Ontology data. This will generate a new folder under `data` containing the necessary data structures for further analysis.
  
```bash
python Data_Processing.py
```

#### 6. Build the Docker Container
- To containerize the application, use the following command to build the Docker image. Replace `<image-name>` with your preferred name for the image.

```bash
docker build -t <image-name> .
```

#### 7. Run the Application
- After building the Docker container, run it using the following command to expose the app on port `8501`:

```bash
docker run -p 8501:8501 <image-name>
```

Once the container is running, you can access the application in your browser at `http://localhost:8501`.

### Notes:
- The AI Assistant uses LangChain agents with built-in tools that map to the three main functionalities.
- The data processing step is crucial to ensure all necessary files and structures are correctly set up for analysis.
- Make sure to have a valid OpenAI API key before running the application.
