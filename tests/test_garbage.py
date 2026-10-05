import errno
import os
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import garbage  # noqa: E402


def make_tree(root: Path, dirs=20, files=5) -> Path:
    for d in range(dirs):
        sub = root / f"d{d}"
        sub.mkdir(parents=True)
        for f in range(files):
            (sub / f"f{f}").write_text("x")
    return root


def wait_empty(root: Path, seconds=20.0) -> bool:
    deadline = time.time() + seconds
    while time.time() < deadline:
        if not [e for e in os.scandir(root) if e.name != garbage.LOCK_NAME]:
            return True
        time.sleep(0.05)
    return False


class GarbageTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.env = mock.patch.dict(os.environ, {}, clear=False)
        self.env.start()
        os.environ.pop(garbage.ENVIRONMENT_VARIABLE, None)

    def tearDown(self):
        self.env.stop()
        import shutil
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_unset_variable_is_plain_rmtree(self):
        work = make_tree(self.tmp / "work")
        with mock.patch.object(garbage.shutil, "rmtree", wraps=garbage.shutil.rmtree) as rmtree:
            garbage.discard(work)
        rmtree.assert_called_once_with(work, ignore_errors=False)
        self.assertFalse(work.exists())

    def test_rename_then_background_delete(self):
        root = garbage.use_garbage(self.tmp / "out" / "garbage")
        work = make_tree(self.tmp / "out" / "work")
        with mock.patch.object(garbage.shutil, "rmtree") as rmtree:
            garbage.discard(work)
        rmtree.assert_not_called()
        self.assertFalse(work.exists())
        self.assertTrue(wait_empty(root), list(os.scandir(root)))

    def test_missing_file_and_symlink_keep_rmtree_behavior(self):
        garbage.use_garbage(self.tmp / "garbage")
        with self.assertRaises(FileNotFoundError):
            garbage.discard(self.tmp / "missing")
        garbage.discard(self.tmp / "missing", ignore_errors=True)
        plain = self.tmp / "file.txt"
        plain.write_text("x")
        with self.assertRaises(NotADirectoryError):
            garbage.discard(plain)
        target = make_tree(self.tmp / "target", dirs=1)
        link = self.tmp / "link"
        link.symlink_to(target)
        with self.assertRaises(OSError):
            garbage.discard(link)
        self.assertTrue(target.exists())

    def test_garbage_inside_discarded_path_deletes_inline(self):
        out = make_tree(self.tmp / "out", dirs=2)
        garbage.use_garbage(out / "garbage")
        garbage.discard(out)
        self.assertFalse(out.exists())

    def test_cross_filesystem_rename_falls_back(self):
        garbage.use_garbage(self.tmp / "garbage")
        work = make_tree(self.tmp / "work", dirs=2)
        with mock.patch.object(garbage.os, "rename", side_effect=OSError(errno.EXDEV, "cross-device")):
            garbage.discard(work)
        self.assertFalse(work.exists())

    def test_many_discards_start_one_deleter(self):
        root = garbage.use_garbage(self.tmp / "garbage")
        works = [make_tree(self.tmp / f"w{i}", dirs=3, files=2) for i in range(30)]
        with mock.patch.object(garbage.subprocess, "Popen") as popen:
            # Simulate a running deleter holding the lock.
            root.mkdir(parents=True, exist_ok=True)
            with open(root / garbage.LOCK_NAME, "a") as handle:
                garbage.fcntl.flock(handle, garbage.fcntl.LOCK_EX)
                for work in works:
                    garbage.discard(work)
            popen.assert_not_called()
        self.assertEqual(len([e for e in os.scandir(root) if e.name != garbage.LOCK_NAME]), 30)
        self.assertEqual(garbage.purge(root), 0)
        self.assertTrue(wait_empty(root))

    def test_sweep_removes_leftovers(self):
        root = self.tmp / "garbage"
        make_tree(root / "left.1.abc", dirs=3)
        garbage.sweep(root)
        self.assertTrue(wait_empty(root))

    def test_purge_stops_when_nothing_can_be_removed(self):
        root = self.tmp / "garbage"
        make_tree(root / "stuck", dirs=1)
        with mock.patch.object(garbage, "_remove_tree"), \
             mock.patch.object(garbage, "start_purger") as restart:
            self.assertEqual(garbage.purge(root), 0)
        restart.assert_not_called()


if __name__ == "__main__":
    unittest.main()
