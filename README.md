## Prerequisites

### Cloud infrastructure (Oracle Cloud)

This project was developed on an Oracle Cloud Infrastructure (OCI) instance.
All pipelines expect a `/data` mount point as the shared working directory.

To set up the data volume:

1. Create a block volume in OCI (paravirtualised) and attach it to your instance
2. Format and mount:
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

3. Make persistent across reboots:
```bash
sudo blkid /dev/sdb   # note the UUID
echo "UUID=<your-uuid> /data ext4 defaults,nofail 0 2" | sudo tee -a /etc/fstab
sudo umount /data
sudo mount -a
df -h /data
```
