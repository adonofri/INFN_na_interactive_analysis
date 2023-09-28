import ROOT

# Apri il file ROOT di input
file_input = ROOT.TFile("/home/jovyan/work/RDataFrame_test/ee_Z_ee_EDM4Hep.root", "READ")

# Controlla se il file e' stato aperto correttamente
if not file_input.IsOpen():
    print("Errore nell'apertura del file di input.")
    exit()

# Leggi il TTree dall'input
tree_input = file_input.Get("events")

# Crea un nuovo file ROOT di output
numero_duplicati = 20  # Specifica quante volte vuoi duplicare il TTree 5,10,20,100
nome_file_output = "input_times_"+str(numero_duplicati)+".root"
file_output = ROOT.TFile(nome_file_output, "RECREATE")
# Crea un nuovo TTree nel file di output con lo stesso nome del TTree originale
tree_duplicato = tree_input.CloneTree(0)

# Ottieni il numero di eventi nel TTree originale
numero_eventi_input = tree_input.GetEntries()

# Duplica il TTree
for _ in range(numero_duplicati):
    tree_duplicato.CopyEntries(tree_input)

# Scrivi il TTree duplicato nel file di output
tree_duplicato.Write()

# Chiudi il file di output
file_output.Close()

# Chiudi il file di input
file_input.Close()
