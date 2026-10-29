# Usa uma imagem oficial do Python baseada em Debian/Ubuntu
FROM python:3.10-slim

# Evita perguntas interativas durante a instalação de pacotes no linux
ENV DEBIAN_FRONTEND=noninteractive

# Instula as ferramentas de sistema necessárias (BLAST, MAFFT, IQ-TREE e compiladores)
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    mafft \
    iqtree \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Define o diretório de trabalho dentro do container
WORKDIR /app

# Copia o arquivo de dependências do Python e instala
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copia todos os arquivos do projeto para dentro do container
COPY . .

# Expõe a porta que o Flask vai rodar
EXPOSE 5000

# Comando para iniciar o servidor usando Gunicorn (recomendado para produção no Render)
CMD gunicorn --bind 0.0.0.0:${PORT:-5000} --timeout 300 --workers 1 servidor_final:app
