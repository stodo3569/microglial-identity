## Prerequisites

### Cloud infrastructure (Oracle Cloud)

This project was developed on Oracle Cloud Infrastructure (OCI) to enable reproducibility, scalability, data security and remote accessibility. The sections below walk through provisioning a virtual machine, configuring networking and security, attaching storage, and deploying the RStudio Server Docker container.

#### 1. Create a virtual machine instance

In the OCI Console, launch a compute instance with the following specifications:

| Setting | Value | Notes |
|---|---|---|
| Shape | VM.Standard.E5.Flex | AMD EPYC; adjustable OCPU and memory |
| OCPU count | 4–8 | 4 minimum for differential expression; 8 recommended for WGCNA |
| Memory (RAM) | 32–80 GB | 32 GB minimum; 64–80 GB recommended for large co-expression analyses |
| Operating system | Ubuntu 24.04 LTS | Canonical image from the OCI marketplace |
| Boot volume | 50 GB (default) | Sufficient for OS + Docker; data goes on a separate block volume |

> **Choosing resources:** The pipeline's most memory-intensive steps are WGCNA network construction and DaMiRseq classification. With fewer than ~15 000 genes these run comfortably in 32 GB, but larger feature sets or multiple concurrent analyses benefit from 64–80 GB. Scale OCPUs and memory to match your dataset size and budget.

During instance creation, upload or paste your **SSH public key**. An SSH key pair consists of two files: a **public key** (installed on the server to verify your identity) and a **private key** (kept securely on your local machine, never shared). OCI needs only the public key; the private key is what you use to authenticate when connecting.

**Connect to the instance via SSH:**

```bash
ssh -i /path/to/private-key ubuntu@<instance-public-ip>
```

The `-i` flag points to your **private key** file. SSH uses it to cryptographically prove your identity to the server, which checks it against the public key you uploaded during instance creation.

#### 2. Configure networking

Create a **Virtual Cloud Network (VCN)** with a public subnet so the instance can reach the internet (required for pulling Docker images and packages):

| Resource | Configuration |
|---|---|
| VCN CIDR | 10.0.0.0/16 (default) |
| Public subnet | 10.0.0.0/24 |
| Internet Gateway | Attached to the VCN |
| Route table rule | Destination 0.0.0.0/0 → Internet Gateway |

#### 3. Configure security

Access control is enforced at two layers — OCI security lists and the OS firewall — to follow a defence-in-depth strategy. OCI security lists act as a virtual firewall at the cloud network level, controlling which traffic is allowed to reach the instance. The OS-level firewall (UFW) provides a second, independent layer of protection running on the instance itself. Even if one layer is misconfigured, the other still blocks unauthorised access.

**OCI Security List (ingress rules):**

"Source" refers to the IP address of the machine initiating the connection. `0.0.0.0/0` means any IP address on the internet. A "port" is a numbered endpoint on the instance that a specific service listens on (e.g. SSH listens on port 22, RStudio Server on port 8787).

| Port | Protocol | Source | Purpose |
|---|---|---|---|
| 22 | TCP | 0.0.0.0/0 | Allow SSH connections from any IP address |
| 8787 | — | *blocked* | RStudio Server — no direct access from the internet; only reachable via an SSH tunnel from within the instance itself (see [step 5](#5-remote-access-to-rstudio-server)) |

**OS-level firewall (UFW):**

```bash
# Allow incoming SSH connections (port 22) so you can log in remotely
sudo ufw allow 22/tcp

# Block any direct connections to RStudio Server (port 8787) from outside;
# it will still be reachable through the SSH tunnel because tunnelled
# traffic arrives as a local (localhost) connection, which UFW does not block
sudo ufw deny 8787/tcp

# Activate the firewall with the rules defined above
sudo ufw enable

# Display the current firewall rules to verify everything is configured correctly
sudo ufw status
```

RStudio Server is **not** exposed to the public internet. Access is restricted to SSH tunnel connections only (see [Remote access to RStudio Server](#5-remote-access-to-rstudio-server) below).

#### 4. Attach and mount a block volume

Create a block volume in OCI (paravirtualised attachment) and attach it to the instance. All pipelines expect `/data` as the shared working directory.

```bash
# Identify the new disk
lsblk

# Format (only once — destroys existing data)
sudo mkfs.ext4 /dev/sdb

# Mount
sudo mkdir -p /data
sudo mount /dev/sdb /data
sudo chown ubuntu:ubuntu /data

# Verify
df -h /data
```

Make persistent across reboots:

```bash
sudo blkid /dev/sdb   # note the UUID
echo "UUID=<your-uuid> /data ext4 defaults,nofail 0 2" | sudo tee -a /etc/fstab
sudo umount /data
sudo mount -a
df -h /data
```

#### 5. Remote access to RStudio Server

RStudio Server listens only on `localhost` inside the VM (see the Docker run command below). To reach it from your local machine, open an SSH tunnel that forwards your local port 8787 to the VM's localhost port 8787:

```bash
# -L 8787:localhost:8787 forwards local port 8787 → VM's localhost:8787
ssh -i /path/to/private-key -L 8787:localhost:8787 ubuntu@<instance-public-ip>
```

Then open **http://localhost:8787** in your browser. The tunnel encrypts all traffic between your machine and the VM.

### Docker image (RStudio Server)

A custom Docker image preserves the complete R / Bioconductor computing environment used for all analyses. The image is based on `bioconductor/bioconductor_docker:RELEASE_3_20` and includes R 4.4.1, RStudio Server, and all required packages (DESeq2, edgeR, sva, DaMiRseq, WGCNA, tidyverse, and others — see `Dockerfile_rstudio-server_microglial-identity` for the full list).

**Pull the image:**

```bash
docker pull stodo3569/rstudio-server_microglial-identity:0.0
```

**Run the container:**

```bash
# -p 127.0.0.1:8787:8787  → Bind RStudio Server to localhost (127.0.0.1) only.
#   This is a universal address meaning "this machine itself". It ensures
#   RStudio is not accessible from the internet; only connections from
#   within the VM (i.e. through the SSH tunnel) can reach port 8787.
# -e DISABLE_AUTH=true     → Skip RStudio login prompt (access is already
#   secured by the SSH tunnel and key-based authentication).
# -v /data:/home/rstudio/data → Mount the host's /data block volume into the
#   container at /home/rstudio/data, so all raw and processed data are
#   accessible within the RStudio environment (host_path:container_path).

docker run -d \
  -p 127.0.0.1:8787:8787 \
  -e DISABLE_AUTH=true \
  -v /data:/home/rstudio/data \
  stodo3569/rstudio-server_microglial-identity:0.0
```

The `-p 127.0.0.1:8787:8787` flag is key to the security model. `127.0.0.1` (localhost) is a standard, universal address — it always refers to the machine itself, regardless of its public IP. By binding to localhost, RStudio Server only accepts connections that originate from within the VM, which is exactly what the SSH tunnel provides.

**Build from source (optional):**

```bash
docker build -t stodo3569/rstudio-server_microglial-identity:0.0 \
  -f Dockerfile_rstudio-server_microglial-identity .
```
