sudo apt-get update
sudo apt install python3.10-venv
sudo apt install openjdk-17-jre -y
echo 'export PATH="/opt/homebrew/opt/
source ~/.zshrc   
java -version
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
nextflow -version



cd /home/ubuntu/data/mydatalocal
git clone https://github.com/Reprohackathon-RNASeq/Reprohackathon-RNASeq.git
cd Reprohackathon-RNASeq/RNASeqanalysis
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt




